// Simplified HB/Elec/Morse pair style
// Author: Margherita
// Modified: HB, Morse, Elec act on ALL non-bonded pairs (intra and inter molecule).
// Intra-molecular HB occupancy is computed geometrically (neighbor scan, same molecule),
// exactly like inter-molecular occupancy. Topology-based occupancy removed.

#include "pair_hb_nonbonded.h"
#include <cmath>
#include <vector>
#include <unordered_map>
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

static const double EPSILON = 1e-12;

/* ---------------------------------------------------------------------- */

PairHBNonbonded::PairHBNonbonded(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  kappa = 0.0;
  dielectric = 1.0;
  use_hb_exclusion = 1;  // Default: exclude electrostatics when HB active
  hb_parallel_threshold = 0.75;   // <-- ADD THIS LINE
}

/* ---------------------------------------------------------------------- */

PairHBNonbonded::~PairHBNonbonded()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(offset);
    memory->destroy(d0);
    memory->destroy(alpha);
    memory->destroy(r0);
    memory->destroy(lambda);
    memory->destroy(epsilon_hb);
    memory->destroy(r_hb);
    memory->destroy(sigma);
    memory->destroy(A_barrier);
  }
}

/* ---------------------------------------------------------------------- */

void PairHBNonbonded::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, ecoul, fpair;
  double rsq, r, rinv, r2inv;
  double u_morse, u_elec, u_hb;
  double f_morse, f_elec;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = ecoul = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  qqrd2e = force->qqrd2e;
  double *special_lj   = force->special_lj;
  double *special_coul = force->special_coul;

  // =====================================================================
  // PRECOMPUTE PER-ATOM HB OCCUPANCY
  // ---------------------------------------------------------------------
  // compute_intra_occ(atom_idx) e compute_inter_occ(atom_idx) dipendono
  // SOLO da atom_idx e dai suoi vicini in neighbor list: NON dipendono
  // dal partner j del pair loop. Prima venivano richiamate una volta per
  // ogni vicino j di i (e viceversa), riscansionando ogni volta l'intera
  // neighbor list di i -> costo O(N*M^2). Qui vengono calcolate una sola
  // volta per atomo e messe in cache -> costo O(N*M).
  // =====================================================================
  int nall = nlocal + atom->nghost;
  std::vector<double> intra_occ_cache(nall, 0.0);
  std::vector<double> inter_occ_cache(nall, 0.0);

  // =====================================================================
  // MAPPA tag -> indice locale (owned + ghost)
  // ---------------------------------------------------------------------
  // is_bonded_12_13_14 doveva prima cercare, per ogni legame 1-3 e 1-4,
  // l'indice locale corrispondente a un tag scorrendo TUTTI gli atomi
  // locali+ghost (O(nall) per ricerca). Con nall grande questa ricerca
  // lineare, ripetuta per ogni vicino di ogni atomo, domina il costo di
  // compute(). Costruendo una volta sola una hash map tag->indice,
  // ogni ricerca diventa O(1) medio invece di O(nall).
  // Costo di costruzione: O(nall), fatto una sola volta per timestep.
  // =====================================================================
  std::unordered_map<tagint,int> tag_to_local;
  tag_to_local.reserve(nall * 2);
  for (int k = 0; k < nall; k++) {
    tag_to_local[atom->tag[k]] = k;
  }

  // Helper: is atom k a 1-2, 1-3 or 1-4 bonded neighbor of atom_idx?
  auto is_bonded_12_13_14 = [&](int atom_idx, int k) -> bool {
    tagint k_tag = atom->tag[k];
    int nb = atom->num_bond[atom_idx];
    for (int ib = 0; ib < nb; ib++) {
      tagint mid_tag = atom->bond_atom[atom_idx][ib];
      if (mid_tag == k_tag) return true;   // 1-2

      // 1-3: check if k is bonded to any neighbor of atom_idx
      auto it_m = tag_to_local.find(mid_tag);
      if (it_m == tag_to_local.end()) continue;  // partner fuori dal dominio locale+ghost
      int m = it_m->second;

      int nb2 = atom->num_bond[m];
      for (int mb = 0; mb < nb2; mb++) {
        tagint mid2_tag = atom->bond_atom[m][mb];
        if (mid2_tag == k_tag) return true;  // 1-3

        // 1-4: check if k is bonded to any 1-3 neighbor
        auto it_p = tag_to_local.find(mid2_tag);
        if (it_p == tag_to_local.end()) continue;
        int p = it_p->second;

        int nb3 = atom->num_bond[p];
        for (int pb = 0; pb < nb3; pb++) {
          if (atom->bond_atom[p][pb] == k_tag) return true;  // 1-4
        }
      }
    }
    return false;
  };

  // ---------------------------------------------------------------
  // GEOMETRY-BASED INTRA-MOLECULAR HB OCCUPANCY
  // Scan neighbor list for atoms of the SAME molecule satisfying
  // HB geometry (distance + normal alignment), exactly as inter.
  // Excludes 1-2 and 1-3 pairs (backbone neighbors).
  // ---------------------------------------------------------------
  auto compute_intra_occ = [&](int atom_idx) -> double {
    double intra_occ = 0.0;
    int total_arr = atom->nlocal + atom->nghost;

    const double d_on  = 7.0;
    const double d_off = 10.0;

    double ni[3];
    compute_normal_vector(x, atom_idx, ni);
    bool has_valid_normal = (atom_idx >= 1 && atom_idx < total_arr - 1);

    if (atom_idx >= list->inum) return 0.0;
    int *jlist_inner = firstneigh[atom_idx];
    int  jnum_inner  = numneigh[atom_idx];

    for (int jj_inner = 0; jj_inner < jnum_inner; jj_inner++) {
      int k = jlist_inner[jj_inner];
      k &= NEIGHMASK;

      // Only atoms of the SAME molecule
      if (atom->molecule[k] != atom->molecule[atom_idx]) continue;

      // Skip 1-2 and 1-3 bonded pairs
      if (is_bonded_12_13_14(atom_idx, k)) continue;

      double dx = x[atom_idx][0] - x[k][0];
      double dy = x[atom_idx][1] - x[k][1];
      double dz = x[atom_idx][2] - x[k][2];
      double dist = sqrt(dx*dx + dy*dy + dz*dz);

      if (dist >= d_off) continue;

      double s_dist;
      if (dist <= d_on) {
        s_dist = 1.0;
      } else {
        double t = (dist - d_on) / (d_off - d_on);
        s_dist = 1.0 - t*t*t*(10.0 - 15.0*t + 6.0*t*t);
      }

      double s_ang = 0.0;
      bool k_interior = (k >= 1 && k < total_arr - 1);
      if (has_valid_normal && k_interior) {
        double nk[3];
        compute_normal_vector(x, k, nk);
        double cos_th = fabs(ni[0]*nk[0] + ni[1]*nk[1] + ni[2]*nk[2]);
        const double ang_low = 0.65, ang_high = 0.75;
        if (cos_th >= ang_high) {
          s_ang = 1.0;
        } else if (cos_th > ang_low) {
          double t = (cos_th - ang_low) / (ang_high - ang_low);
          s_ang = t*t*(3.0 - 2.0*t);
        }
      }
      intra_occ += s_dist * s_ang;
    }
    if (intra_occ > 2.0) intra_occ = 2.0;
    return intra_occ;
  };

  // ---------------------------------------------------------------
  // GEOMETRY-BASED INTER-MOLECULAR HB OCCUPANCY (unchanged logic)
  // Scan neighbor list for atoms of a DIFFERENT molecule.
  // ---------------------------------------------------------------
  auto compute_inter_occ = [&](int atom_idx) -> double {
    double inter_occ = 0.0;
    int total_arr = atom->nlocal + atom->nghost;

    const double d_on  = 7.0;
    const double d_off = 10.0;

    double ni[3];
    compute_normal_vector(x, atom_idx, ni);
    bool has_valid_normal = (atom_idx >= 1 && atom_idx < total_arr - 1);

    if (atom_idx >= list->inum) return 0.0;
    int *jlist_inner = firstneigh[atom_idx];
    int  jnum_inner  = numneigh[atom_idx];

    for (int jj_inner = 0; jj_inner < jnum_inner; jj_inner++) {
      int k = jlist_inner[jj_inner];
      k &= NEIGHMASK;

      // Only atoms of a DIFFERENT molecule
      if (atom->molecule[k] == atom->molecule[atom_idx]) continue;

      double dx = x[atom_idx][0] - x[k][0];
      double dy = x[atom_idx][1] - x[k][1];
      double dz = x[atom_idx][2] - x[k][2];
      double dist = sqrt(dx*dx + dy*dy + dz*dz);

      if (dist >= d_off) continue;

      double s_dist;
      if (dist <= d_on) {
        s_dist = 1.0;
      } else {
        double t = (dist - d_on) / (d_off - d_on);
        s_dist = 1.0 - t*t*t*(10.0 - 15.0*t + 6.0*t*t);
      }

      double s_ang = 0.0;
      bool k_interior = (k >= 1 && k < total_arr - 1);
      if (has_valid_normal && k_interior) {
        double nk[3];
        compute_normal_vector(x, k, nk);
        double cos_th = fabs(ni[0]*nk[0] + ni[1]*nk[1] + ni[2]*nk[2]);
        const double ang_low = 0.65, ang_high = 0.75;
        if (cos_th >= ang_high) {
          s_ang = 1.0;
        } else if (cos_th > ang_low) {
          double t = (cos_th - ang_low) / (ang_high - ang_low);
          s_ang = t*t*(3.0 - 2.0*t);
        }
      }
      inter_occ += s_dist * s_ang;
    }
    if (inter_occ > 2.0) inter_occ = 2.0;
    return inter_occ;
  };

  // Riempi la cache UNA volta per ciascun atomo owned (i in ilist).
  // Gli atomi ghost restano a 0.0, esattamente come faceva il codice
  // originale (return 0.0 per atom_idx >= list->inum).
  for (int ii_pre = 0; ii_pre < inum; ii_pre++) {
    int iatom = ilist[ii_pre];
    intra_occ_cache[iatom] = compute_intra_occ(iatom);
    inter_occ_cache[iatom] = compute_inter_occ(iatom);
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      double factor_lj   = special_lj[sbmask(j)];
      double factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        rinv = 1.0 / r;
        r2inv = rinv * rinv;

        fpair = 0.0;
        u_morse = 0.0;
        u_elec = 0.0;
        u_hb = 0.0;
        f_morse = 0.0;
        f_elec = 0.0;
        
        // HB force vector (non-radial)
        double f_hb[3] = {0.0, 0.0, 0.0};

        // ========== MORSE POTENTIAL (RADIAL) ==========
        if (d0[itype][jtype] > 0.0) {
          double dr = r - r0[itype][jtype];
          double dexp = exp(-alpha[itype][jtype] * dr);
          double dexp2 = dexp * dexp;
          
          u_morse = d0[itype][jtype] * (dexp2 - 2.0 * dexp);
          f_morse = 2.0 * d0[itype][jtype] * alpha[itype][jtype] * (dexp2 - dexp) * rinv;
          u_morse *= factor_lj;
          f_morse *= factor_lj;

        }

        // ========== ELECTROSTATIC POTENTIAL (RADIAL) ==========
        if (lambda[itype][jtype] != 0.0) {
          double screening = (kappa > 0.0) ? exp(-kappa * r) : 1.0;
          u_elec = qqrd2e * lambda[itype][jtype] * screening / (dielectric * r);
          double f_elec_factor = screening * (kappa + rinv) * r2inv;
          f_elec = qqrd2e * lambda[itype][jtype] * f_elec_factor / dielectric;
          u_elec *= factor_coul;
          f_elec *= factor_coul;

        }



        // ========== HYDROGEN BOND POTENTIAL (NON-RADIAL) ==========
        if (epsilon_hb[itype][jtype] != 0.0) {

          bool same_molecule = (atom->molecule[i] == atom->molecule[j]);

          // Returns true if i and j are connected by a declared HB bond (bond type > 1).
          // These pairs are excluded from the HB potential (no double-counting with bonded),
          // but ARE included in compute_intra_occ to contribute to occupancy.
          auto is_hb_bond = [&](int atom_idx, int k) -> bool {
            tagint k_tag = atom->tag[k];
            int nb = atom->num_bond[atom_idx];
            int *btype = atom->bond_type[atom_idx];
            tagint *batom = atom->bond_atom[atom_idx];
            for (int ib = 0; ib < nb; ib++) {
              if (btype[ib] > 1 && batom[ib] == k_tag) return true;
            }
            return false;
          };

          // ---------------------------------------------------------------
          // Lookup dell'occupancy dalla cache pre-calcolata (una volta per
          // atomo, PRIMA del doppio loop sulle coppie). Nessun ricalcolo
          // qui: i valori sono già pronti in intra_occ_cache/inter_occ_cache.
          // Per gli atomi ghost (j >= inum) il valore resta 0.0, identico
          // al comportamento originale.
          // ---------------------------------------------------------------
          double i_hb_occupancy = intra_occ_cache[i];
          double j_hb_occupancy = intra_occ_cache[j];
          double i_inter_occ    = inter_occ_cache[i];
          double j_inter_occ    = inter_occ_cache[j];

          // ---------------------------------------------------------------
          // MUTUAL EXCLUSION:
          // (A) Inter suppresses intra: reduce intra occupancy by inter occupancy.
          //     When an inter-HB is formed, the intra occupancy is partially reduced
          //     so it does not block the inter HB.
          // (B) Intra suppresses inter: when intra-HB occupancy is high, the HB
          //     potential for the current pair (inter or intra) is scaled down.
          //
          // For intra pairs: inter_hb_scale is always 1.0 (intra does not suppress
          // itself via the inter gate — the intra occ of i and j is the same pool).
          // ---------------------------------------------------------------
          auto ramp1 = [](double occ) -> double {
            if (occ <= 0.0) return 0.0;
            if (occ >= 1.0) return 1.0;
            double u = occ;
            return u * u * (3.0 - 2.0 * u);
          };

          // (A) Reduce intra occupancy by intermolecular occupancy.
          const double alpha_suppress = 1.00;
          double inter_suppress_factor_i = 1.0 - alpha_suppress * ramp1(i_inter_occ);
          double inter_suppress_factor_j = 1.0 - alpha_suppress * ramp1(j_inter_occ);
          i_hb_occupancy *= inter_suppress_factor_i;
          j_hb_occupancy *= inter_suppress_factor_j;

          // (B) Smooth gate: suppress inter-HB when intra occupancy is high.
          // For same-molecule pairs, no self-suppression via this gate.
          //
          // Tuning vs original:
          //   threshold raised 0.5 -> 0.8: gate only activates when intra nearly saturated
          //   max suppression reduced: 1.0 -> 0.5: inter HB never fully blocked by intra
          auto ramp2 = [](double occ) -> double {
            const double low = 0.8;   // was 0.5
            if (occ <= low) return 0.0;
            if (occ >= 1.0) return 1.0;
            double u = (occ - low) / (1.0 - low);
            return u * u * (3.0 - 2.0 * u);
          };
          double inter_hb_scale;
          if (same_molecule) {
            // Intra pair: no inter-gate suppression; intra occ does not suppress itself here.
            inter_hb_scale = 1.0;
          } else {
            double suppress_i = ramp2(i_hb_occupancy);
            double suppress_j = ramp2(j_hb_occupancy);
            double max_suppress = (suppress_i > suppress_j) ? suppress_i : suppress_j;
            inter_hb_scale = 1.0 - 0.5 * max_suppress;  // was 1.0*max_suppress
          }

          // ---------------------------------------------------------------
          // ---------------------------------------------------------------
          // REPULSIVE PENALTY: when inter-HB is formed (intra_suppress > 0.5),
          // apply a soft repulsion between i and j (the current intra pair)
          // to destabilize the intra contact. This acts directly on the pair
          // being evaluated, not on pre-declared bond partners.
          // Only applied to same-molecule pairs.
          // ---------------------------------------------------------------
          double intra_suppress_i = ramp1(i_inter_occ);
          double intra_suppress_j = ramp1(j_inter_occ);

          if (same_molecule) {
            double intra_suppress = (intra_suppress_i + intra_suppress_j) * 0.5;
            if (intra_suppress > 0.5) {
              const double eps_rep = 0.07, r_rep = 6.5, sig_rep = 0.8;
              double dr_rep = r - r_rep;
              double sigma_sq = sig_rep * sig_rep;
              double u_rep = eps_rep * intra_suppress * exp(-dr_rep*dr_rep / sigma_sq);
              double f_rep_mag = u_rep * (2.0 * dr_rep / sigma_sq) * rinv;
              // accumulate into fpair (radial, same sign convention as Morse)
              fpair += f_rep_mag;
              if (eflag) evdwl += u_rep;
            }
          }

          // ---------------------------------------------------------------
          // HB POTENTIAL: active for pairs that are NOT:
          //   - 1-2 or 1-3 bonded (backbone neighbors), AND
          //   - already connected by a declared HB bond (bond type > 1).
          // HB bonds are excluded here to avoid double-counting with the
          // bonded part, but they DO contribute to compute_intra_occ above.
          // Same molecule or different molecule — geometry decides.
          // ---------------------------------------------------------------
          bool are_bonded_12_13_14 = is_bonded_12_13_14(i, j);
          bool are_hb_bonded    = is_hb_bond(i, j);
          if (!are_bonded_12_13_14 && !are_hb_bonded) {
            
            // Guard: compute_normal_vector needs i-1, i, i+1 and j-1, j, j+1 to be valid.
            // Use total array size (local + ghost) as the bound.
            int total_with_ghosts = nlocal + atom->nghost;
            bool can_compute_hb = true;
            if (i < 1 || i >= total_with_ghosts - 1 || j < 1 || j >= total_with_ghosts - 1) {
              can_compute_hb = false;
            }
            
            if (can_compute_hb) {
              double rij_hat[3];
              rij_hat[0] = delx * rinv;
              rij_hat[1] = dely * rinv;
              rij_hat[2] = delz * rinv;

              double ti_unit[3], tj_unit[3];
              compute_normal_vector(x, i, ti_unit);
              compute_normal_vector(x, j, tj_unit);

              double alpha_i = rij_hat[0] * ti_unit[0] + rij_hat[1] * ti_unit[1] + rij_hat[2] * ti_unit[2];
              double alpha_j = rij_hat[0] * tj_unit[0] + rij_hat[1] * tj_unit[1] + rij_hat[2] * tj_unit[2];

              // Use unified sigma for both distance and angular
              double sigma_sq = sigma[itype][jtype] * sigma[itype][jtype];
              
              // Energy components
              double dr_hb = r - r_hb[itype][jtype];
              double h0 = -epsilon_hb[itype][jtype] * exp(-dr_hb * dr_hb / sigma_sq);
              double h1 = exp((fabs(alpha_i) - 1.0) / sigma_sq);
              double h2 = exp((fabs(alpha_j) - 1.0) / sigma_sq);
              
              // ---- Soft-capped repulsive barrier ----
              // The raw barrier ~ (r_hb/r)^12 diverges as r -> 0 (force ~1/r^13,
              // steeper than Morse or electrostatics). A transient close contact
              // during an HB forming/breaking event can push r into this regime
              // faster than the timestep can resolve, giving a force spike large
              // enough to move an atom past comm_modify cutoff in one step -> the
              // "Bond atoms missing" crash. Below r_cap we switch to a linear
              // (constant-force) extrapolation that matches value and slope at
              // r_cap, so the potential stays C1-continuous and the force is
              // bounded. r_cap_ratio=0.75 leaves the wall untouched in the
              // physically relevant range and only clips the pathological regime;
              // tune if it turns out to bite during normal HB formation.
              const double r_cap_ratio = 0.75;
              double r_cap = r_cap_ratio * r_hb[itype][jtype];

              double barrier, dbarrier_dr_capped;
              if (r >= r_cap) {
                double r_shifted = r * pow(2.0, 1.0/6.0);  // r_effective = 2^(1/6) * r
                double r_ratio = r_hb[itype][jtype] / r_shifted;
                barrier = A_barrier[itype][jtype] * pow(r_ratio, 12.0);
                dbarrier_dr_capped = -12.0 * barrier * rinv;
              } else {
                double r_shifted_cap = r_cap * pow(2.0, 1.0/6.0);
                double r_ratio_cap = r_hb[itype][jtype] / r_shifted_cap;
                double barrier_cap = A_barrier[itype][jtype] * pow(r_ratio_cap, 12.0);
                double dbarrier_dr_cap = -12.0 * barrier_cap / r_cap;  // constant (bounded) force below r_cap
                barrier = barrier_cap + dbarrier_dr_cap * (r - r_cap);
                dbarrier_dr_capped = dbarrier_dr_cap;
              }
              
              u_hb = (h0 * h1 * h2 + barrier) * inter_hb_scale * factor_lj;

              // ========== HB-ELECTROSTATICS EXCLUSION ==========
              if (use_hb_exclusion && fabs(lambda[itype][jtype]) > EPSILON) {
                // Normalize HB strength (0 to 1): when HB is strong, turn off electrostatics
                double hb_normalized = fabs(h0 * h1 * h2) / epsilon_hb[itype][jtype];
                if (hb_normalized > 1.0) hb_normalized = 1.0;
                
                // Quintic switching function: smooth transition from 1 (no HB) to 0 (full HB)
                double x = hb_normalized;
                double switch_elec = 1.0 - x * x * x * (10.0 - 15.0 * x + 6.0 * x * x);
                
                u_elec *= switch_elec;
                f_elec *= switch_elec;
              }

              // ========== FORCE CALCULATION (NON-RADIAL) ==========
              double dh0_dr = h0 * (-2.0 * dr_hb / sigma_sq);
              double dbarrier_dr = dbarrier_dr_capped;  // bounded, see capping block above
              
              double sgn_alpha_i = (fabs(alpha_i) < EPSILON) ? 0.0 : ((alpha_i > 0.0) ? 1.0 : -1.0);
              double sgn_alpha_j = (fabs(alpha_j) < EPSILON) ? 0.0 : ((alpha_j > 0.0) ? 1.0 : -1.0);
              
              // Angular derivatives (non-radial components)
              double dh1_coeff = (h1 * sgn_alpha_i) / (sigma_sq * r);
              double dh1_dri[3];
              for (int d = 0; d < 3; d++) {
                dh1_dri[d] = dh1_coeff * (ti_unit[d] - alpha_i * rij_hat[d]);
              }
              
              double dh2_coeff = (h2 * sgn_alpha_j) / (sigma_sq * r);
              double dh2_dri[3];
              for (int d = 0; d < 3; d++) {
                dh2_dri[d] = dh2_coeff * (tj_unit[d] - alpha_j * rij_hat[d]);
              }
              
              // Total HB force: F = -dU/dr, scaled by inter_hb_scale
              for (int d = 0; d < 3; d++) {
                    f_hb[d] = -((dh0_dr * rij_hat[d]) * h1 * h2
                              + h0 * dh1_dri[d] * h2
                              + h0 * h1 * dh2_dri[d]
                              + dbarrier_dr * rij_hat[d])
                              * inter_hb_scale * factor_lj;
                }                             // smooth suppression
            }  // end can_compute_hb
          }  // end !are_bonded_12_13 && !are_hb_bonded
        }  // end epsilon_hb check

        // ========== APPLY FORCES ==========
        // Radial forces (Morse + Electrostatics + repulsive penalty): F_radial = fpair * r_hat
        fpair += f_morse + f_elec;
        
        // Apply radial forces
        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        
        // Apply non-radial HB forces
        f[i][0] += f_hb[0];
        f[i][1] += f_hb[1];
        f[i][2] += f_hb[2];
        
        if (newton_pair || j < nlocal) {
          // Newton's third law
          f[j][0] -= delx * fpair + f_hb[0];
          f[j][1] -= dely * fpair + f_hb[1];
          f[j][2] -= delz * fpair + f_hb[2];
        }

        if (eflag) {
          evdwl = u_morse + u_hb - offset[itype][jtype];
          ecoul = u_elec;
        }

        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, ecoul, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairHBNonbonded::compute_normal_vector(double **x, int i, double *t_unit)
{
  // Fallback vector
  auto fallback = [&]() {
    t_unit[0] = 0.0;
    t_unit[1] = 0.0;
    t_unit[2] = 1.0;
  };

  int total = atom->nlocal + atom->nghost;
  tagint *tag = atom->tag;

  // Collect indices of ALL bonded neighbors (harmonic bonds, type 1 = backbone)
  // We want exactly the two backbone neighbors to define the local plane
  int prev_idx = -1, next_idx = -1;
  int num_b = atom->num_bond[i];

  for (int ib = 0; ib < num_b; ib++) {
    // Only use backbone (harmonic) bonds, type == 1
    if (atom->bond_type[i][ib] != 1) continue;

    tagint bonded_tag = atom->bond_atom[i][ib];

    // Search for this tag in local + ghost array
    for (int k = 0; k < total; k++) {
      if (tag[k] == bonded_tag) {
        if (prev_idx < 0) prev_idx = k;
        else              next_idx = k;
        break;
      }
    }
  }

  // Need exactly 2 backbone neighbors to define a plane
  if (prev_idx < 0 || next_idx < 0) {
    fallback();
    return;
  }

  // Cross product of (i -> prev) x (i -> next)
  double vi1[3], vi2[3], ti[3];
  for (int d = 0; d < 3; d++) {
    vi1[d] = x[i][d] - x[prev_idx][d];
    vi2[d] = x[next_idx][d] - x[i][d];
  }

  ti[0] = vi1[1]*vi2[2] - vi1[2]*vi2[1];
  ti[1] = vi1[2]*vi2[0] - vi1[0]*vi2[2];
  ti[2] = vi1[0]*vi2[1] - vi1[1]*vi2[0];

  double ti_norm = sqrt(ti[0]*ti[0] + ti[1]*ti[1] + ti[2]*ti[2]);

  if (ti_norm > EPSILON) {
    t_unit[0] = ti[0] / ti_norm;
    t_unit[1] = ti[1] / ti_norm;
    t_unit[2] = ti[2] / ti_norm;
  } else {
    fallback();
  }
}

/* ---------------------------------------------------------------------- */

void PairHBNonbonded::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(offset, n + 1, n + 1, "pair:offset");
  memory->create(d0, n + 1, n + 1, "pair:d0");
  memory->create(alpha, n + 1, n + 1, "pair:alpha");
  memory->create(r0, n + 1, n + 1, "pair:r0");
  memory->create(lambda, n + 1, n + 1, "pair:lambda");
  memory->create(epsilon_hb, n + 1, n + 1, "pair:epsilon_hb");
  memory->create(r_hb, n + 1, n + 1, "pair:r_hb");
  memory->create(sigma, n + 1, n + 1, "pair:sigma");
  memory->create(A_barrier, n + 1, n + 1, "pair:A_barrier");

  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;
}

/* ---------------------------------------------------------------------- */

void PairHBNonbonded::settings(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR, "Illegal pair_style command");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);
  if (narg >= 2) kappa = utils::numeric(FLERR, arg[1], false, lmp);
  if (narg >= 3) dielectric = utils::numeric(FLERR, arg[2], false, lmp);
  
  if (allocated) {
    int n = atom->ntypes;
    for (int i = 1; i <= n; i++) {
      for (int j = i; j <= n; j++) {
        if (setflag[i][j] == 0) cut[i][j] = cut_global;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Simplified coeff format:
   pair_coeff i j D0 alpha r0 lambda epsilon_hb r_hb sigma A [cutoff]
------------------------------------------------------------------------- */

void PairHBNonbonded::coeff(int narg, char **arg)
{
  if (narg < 10 && narg > 11)
    error->all(FLERR, "Incorrect args for pair coefficients (need 8 params + optional cutoff)");
  
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double d0_one = utils::numeric(FLERR, arg[2], false, lmp);
  double alpha_one = utils::numeric(FLERR, arg[3], false, lmp);
  double r0_one = utils::numeric(FLERR, arg[4], false, lmp);
  double lambda_one = utils::numeric(FLERR, arg[5], false, lmp);
  double epsilon_hb_one = utils::numeric(FLERR, arg[6], false, lmp);
  double r_hb_one = utils::numeric(FLERR, arg[7], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[8], false, lmp);
  double A_one = utils::numeric(FLERR, arg[9], false, lmp);

  double cut_one = cut_global;
  if (narg == 11) cut_one = utils::numeric(FLERR, arg[10], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = jlo; j <= jhi; j++) {
      d0[i][j] = d0_one;
      alpha[i][j] = alpha_one;
      r0[i][j] = r0_one;
      lambda[i][j] = lambda_one;
      epsilon_hb[i][j] = epsilon_hb_one;
      r_hb[i][j] = r_hb_one;
      sigma[i][j] = sigma_one;
      A_barrier[i][j] = A_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ---------------------------------------------------------------------- */

void PairHBNonbonded::init_style()
{
  if (!atom->molecule_flag)
    error->all(FLERR, "Pair style hb/simple requires molecule IDs");
  neighbor->request(this, instance_me);
}

/* ---------------------------------------------------------------------- */

double PairHBNonbonded::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  offset[i][j] = 0.0;
  if (offset_flag) {
    double r = cut[i][j];
    double rinv = 1.0 / r;
    
    if (d0[i][j] > 0.0) {
      double dr = r - r0[i][j];
      double dexp = exp(-alpha[i][j] * dr);
      offset[i][j] += d0[i][j] * (dexp * dexp - 2.0 * dexp);
    }
    
    if (epsilon_hb[i][j] > 0.0) {
      double sigma_sq = sigma[i][j] * sigma[i][j];
      double dr_hb = r - r_hb[i][j];
      double h0 = -epsilon_hb[i][j] * exp(-dr_hb * dr_hb / sigma_sq);
      double r_shifted = r * pow(2.0, 1.0/6.0);
      double barrier = A_barrier[i][j] * pow(r_hb[i][j] / r_shifted, 12.0);
      offset[i][j] += h0 + barrier;
    }
  }

  d0[j][i] = d0[i][j];
  alpha[j][i] = alpha[i][j];
  r0[j][i] = r0[i][j];
  lambda[j][i] = lambda[i][j];
  epsilon_hb[j][i] = epsilon_hb[i][j];
  r_hb[j][i] = r_hb[i][j];
  sigma[j][i] = sigma[i][j];
  A_barrier[j][i] = A_barrier[i][j];
  cut[j][i] = cut[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

void PairHBNonbonded::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&d0[i][j], sizeof(double), 1, fp);
        fwrite(&alpha[i][j], sizeof(double), 1, fp);
        fwrite(&r0[i][j], sizeof(double), 1, fp);
        fwrite(&lambda[i][j], sizeof(double), 1, fp);
        fwrite(&epsilon_hb[i][j], sizeof(double), 1, fp);
        fwrite(&r_hb[i][j], sizeof(double), 1, fp);
        fwrite(&sigma[i][j], sizeof(double), 1, fp);
        fwrite(&A_barrier[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairHBNonbonded::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &d0[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &alpha[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &r0[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &lambda[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &epsilon_hb[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &r_hb[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &sigma[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &A_barrier[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&d0[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&alpha[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&r0[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&lambda[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&epsilon_hb[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&r_hb[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&A_barrier[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairHBNonbonded::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&kappa, sizeof(double), 1, fp);
  fwrite(&dielectric, sizeof(double), 1, fp);
  fwrite(&use_hb_exclusion, sizeof(int), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ---------------------------------------------------------------------- */

void PairHBNonbonded::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &kappa, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &dielectric, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &use_hb_exclusion, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&kappa, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dielectric, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&use_hb_exclusion, 1, MPI_INT, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

void PairHBNonbonded::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g %g %g %g %g %g %g %g\n", i,
            d0[i][i], alpha[i][i], r0[i][i], lambda[i][i],
            epsilon_hb[i][i], r_hb[i][i], sigma[i][i], A_barrier[i][i]);
}

/* ---------------------------------------------------------------------- */

void PairHBNonbonded::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      fprintf(fp, "%d %d %g %g %g %g %g %g %g %g %g\n", i, j,
              d0[i][j], alpha[i][j], r0[i][j], lambda[i][j],
              epsilon_hb[i][j], r_hb[i][j], sigma[i][j], A_barrier[i][j], cut[i][j]);
    }
  }
}
