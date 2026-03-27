#include "pair_hb_nonbonded.h"
#include <cmath>
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
          
          // Check if atoms are from different molecules using LAMMPS molecule IDs.
          // Works for any number of chains, as long as mol-id is set correctly
          // in the data file Atoms section.
          tagint mol_i = atom->molecule[i];
          tagint mol_j = atom->molecule[j];
          bool different_molecules = (mol_i != mol_j);
          
          // SMOOTH VERSION: compute continuous occupancy of intramolecular hb6barriernew bonds
          // instead of a binary count. occupancy in [0,1] per bond; sum capped at 1.
          // This avoids abrupt force discontinuities when the intra-HB breaks/forms.
          double i_hb_occupancy = 0.0;
          int num_bonds_i = atom->num_bond[i];
          if (num_bonds_i > 0) {
            int *bond_type_i = atom->bond_type[i];
            int *bond_atom_i = atom->bond_atom[i];
            tagint *tag = atom->tag;
            int nlocal_atoms = atom->nlocal;
            int total_atoms_arr = nlocal_atoms + atom->nghost;

            // Distance switching: 1 at r=0, smoothly 0 at r=d_switch_off
            const double d_switch_on  = 7.0;   // fully formed below this distance
            const double d_switch_off = 10.0;   // fully broken above this distance

            double ni[3];
            compute_normal_vector(x, i, ni);
            bool i_has_valid_normal = (i >= 1 && i < total_atoms_arr - 1);
            
            for (int ib = 0; ib < num_bonds_i; ib++) {
              if (bond_type_i[ib] > 1) {
                tagint bonded_tag = bond_atom_i[ib];
                
                int bonded_idx = -1;
                for (int k = 0; k < total_atoms_arr; k++) {
                  if (tag[k] == bonded_tag) { bonded_idx = k; break; }
                }
                
                if (bonded_idx >= 0 && atom->molecule[bonded_idx] == atom->molecule[i]) {
                  double dx = x[i][0] - x[bonded_idx][0];
                  double dy = x[i][1] - x[bonded_idx][1];
                  double dz = x[i][2] - x[bonded_idx][2];
                  double bond_dist = sqrt(dx*dx + dy*dy + dz*dz);
                  
                  if (bond_dist < d_switch_off) {
                    // Smooth 5th-order distance switch: 1 when close, 0 when far
                    double s_dist;
                    if (bond_dist <= d_switch_on) {
                      s_dist = 1.0;
                    } else {
                      double t = (bond_dist - d_switch_on) / (d_switch_off - d_switch_on); // 0->1
                      s_dist = 1.0 - t * t * t * (10.0 - 15.0 * t + 6.0 * t * t);
                    }
                    
                    // Smooth angular switch: 1 when normals are aligned, 0 at threshold
                    double s_ang = 0.0;
                    bool bonded_is_interior = (bonded_idx >= 1 && bonded_idx < total_atoms_arr - 1);
                    if (i_has_valid_normal && bonded_is_interior) {
                        double nb[3];
                        compute_normal_vector(x, bonded_idx, nb);
                        double cos_theta = fabs(ni[0]*nb[0] + ni[1]*nb[1] + ni[2]*nb[2]);
                        const double ang_low  = 0.65;
                        const double ang_high = 0.75;  // == hb_parallel_threshold
                        if (cos_theta >= ang_high) {
                          s_ang = 1.0;
                        } else if (cos_theta > ang_low) {
                          double t = (cos_theta - ang_low) / (ang_high - ang_low);
                          s_ang = t * t * (3.0 - 2.0 * t);
                        }
                      }
                    
                    i_hb_occupancy += s_dist * s_ang;
                  }
                }
              }
            }
            if (i_hb_occupancy > 2.0) i_hb_occupancy = 2.0;  // cap at 2: max 2 intra-HBs counted
          }
          
          // SMOOTH VERSION: same continuous occupancy for atom j
          double j_hb_occupancy = 0.0;
          int num_bonds_j = atom->num_bond[j];
          if (num_bonds_j > 0) {
            int *bond_type_j = atom->bond_type[j];
            int *bond_atom_j = atom->bond_atom[j];
            tagint *tag = atom->tag;
            int nlocal_atoms = atom->nlocal;
            int total_atoms_arr = nlocal_atoms + atom->nghost;

            const double d_switch_on  = 7.0;
            const double d_switch_off = 10.0;

            double nj[3];
            compute_normal_vector(x, j, nj);
            bool j_has_valid_normal = (j >= 1 && j < total_atoms_arr - 1);
            
            for (int jb = 0; jb < num_bonds_j; jb++) {
              if (bond_type_j[jb] > 1) {
                tagint bonded_tag = bond_atom_j[jb];
                
                int bonded_idx = -1;
                for (int k = 0; k < total_atoms_arr; k++) {
                  if (tag[k] == bonded_tag) { bonded_idx = k; break; }
                }
                
                if (bonded_idx >= 0 && atom->molecule[bonded_idx] == atom->molecule[j]) {
                  double dx = x[j][0] - x[bonded_idx][0];
                  double dy = x[j][1] - x[bonded_idx][1];
                  double dz = x[j][2] - x[bonded_idx][2];
                  double bond_dist = sqrt(dx*dx + dy*dy + dz*dz);
                  
                  if (bond_dist < d_switch_off) {
                    double s_dist;
                    if (bond_dist <= d_switch_on) {
                      s_dist = 1.0;
                    } else {
                      double t = (bond_dist - d_switch_on) / (d_switch_off - d_switch_on);
                      s_dist = 1.0 - t * t * t * (10.0 - 15.0 * t + 6.0 * t * t);
                    }
                    
                    double s_ang = 0.0;
                    bool bonded_is_interior = (bonded_idx >= 1 && bonded_idx < total_atoms_arr - 1);
                    if (j_has_valid_normal && bonded_is_interior) {
                        double nb[3];
                        compute_normal_vector(x, bonded_idx, nb);
                        double cos_theta = fabs(nj[0]*nb[0] + nj[1]*nb[1] + nj[2]*nb[2]);
                        const double ang_low  = 0.65;
                        const double ang_high = 0.75;
                        if (cos_theta >= ang_high) {
                          s_ang = 1.0;
                        } else if (cos_theta > ang_low) {
                          double t = (cos_theta - ang_low) / (ang_high - ang_low);
                          s_ang = t * t * (3.0 - 2.0 * t);
                        }
                      }
                    
                    j_hb_occupancy += s_dist * s_ang;
                  }
                }
              }
            }
            if (j_hb_occupancy > 2.0) j_hb_occupancy = 2.0;  // cap at 2: max 2 intra-HBs counted
          }

          // ---------------------------------------------------------------
          // COMPUTE INTERMOLECULAR HB OCCUPANCY for atoms i and j
          // Scans ALL atoms of the other molecule in the neighbor list
          // (geometry-based, no pre-declared bond required)
          // ---------------------------------------------------------------
          auto compute_inter_occ = [&](int atom_idx, tagint mol_id) -> double {
            double inter_occ = 0.0;
            int total_arr = atom->nlocal + atom->nghost;

            const double d_on  = 7.0;
            const double d_off = 10.0;

            double ni[3];
            compute_normal_vector(x, atom_idx, ni);
            bool has_valid_normal = (atom_idx >= 1 && atom_idx < total_arr - 1);

            // Use neighbor list: scan all neighbors of atom_idx
            // firstneigh/numneigh are populated for all atoms in ilist;
            // for ghost atoms (j) the list may be empty, which is safe.
            if (atom_idx >= list->inum) return 0.0;  // atom_idx not in ilist, skip
            int *jlist_inner = firstneigh[atom_idx];
            int  jnum_inner  = numneigh[atom_idx];

            for (int jj_inner = 0; jj_inner < jnum_inner; jj_inner++) {
              int k = jlist_inner[jj_inner];
              k &= NEIGHMASK;

              // Only consider atoms from a DIFFERENT molecule (works for any number of chains)
              if (atom->molecule[k] == mol_id) continue;  // same molecule, skip

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
                s_dist = 1.0 - t * t * t * (10.0 - 15.0 * t + 6.0 * t * t);
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
                  s_ang = t * t * (3.0 - 2.0 * t);
                }
              }
              inter_occ += s_dist * s_ang;
            }
            if (inter_occ > 2.0) inter_occ = 2.0;
            return inter_occ;
          };

          double i_inter_occ = compute_inter_occ(i, mol_i);
          double j_inter_occ = compute_inter_occ(j, mol_j);

          // ---------------------------------------------------------------
          // MUTUAL EXCLUSION:
          // (A) Inter suppresses intra: reduce intra occupancy by inter occupancy
          //     so that when inter is formed, intra is not counted as blocking inter.
          // (B) Intra suppresses inter: original logic below.
          // ---------------------------------------------------------------
          // Smooth ramp: 0 at occ=0, 1 at occ>=1
          auto ramp1 = [](double occ) -> double {
            if (occ <= 0.0) return 0.0;
            if (occ >= 1.0) return 1.0;
            double u = occ;
            return u * u * (3.0 - 2.0 * u);
          };

          // (A) Reduce intra occupancy by intermolecular occupancy (mutual exclusion)
          double inter_suppress_factor_i = 1.0 - ramp1(i_inter_occ);
          double inter_suppress_factor_j = 1.0 - ramp1(j_inter_occ);
          i_hb_occupancy *= inter_suppress_factor_i;
          j_hb_occupancy *= inter_suppress_factor_j;

          // ---------------------------------------------------------------
          // SMOOTH GATE: suppress intermolecular HB when intra occupancy is high.
          // (original logic, now with intra occ already reduced by inter occ above)
          // ---------------------------------------------------------------
          auto ramp2 = [](double occ) -> double {
            if (occ <= 0.5) return 0.0;
            if (occ >= 1.0) return 1.0;
            double u = (occ - 0.5) / 0.5;
            return u * u * (3.0 - 2.0 * u);
          };
          double suppress_i = ramp2(i_hb_occupancy);
          double suppress_j = ramp2(j_hb_occupancy);
          double max_suppress = (suppress_i > suppress_j) ? suppress_i : suppress_j;
          double inter_hb_scale = 1.0 - max_suppress;

          // ---------------------------------------------------------------
          // (B) INTRA SUPPRESSION BY INTER: add repulsive penalty on the
          //     intramolecular HB bonded pairs of i and j when inter occ is high.
          //     This directly destabilizes the intra HB when inter is formed.
          //     We apply it as an additional energy/force contribution here.
          // ---------------------------------------------------------------
          // Scale of repulsion: when inter_occ >= 1, apply full repulsion
          double intra_suppress_i = ramp1(i_inter_occ);  // 0 when no inter, 1 when inter formed
          double intra_suppress_j = ramp1(j_inter_occ);

          // Apply repulsive penalty to intramolecular HB bond partners of i
          if (intra_suppress_i > 1e-6) {
            int num_b_i = atom->num_bond[i];
            int *btype_i = atom->bond_type[i];
            int *batom_i = atom->bond_atom[i];
            tagint *tag = atom->tag;
            int total_arr = atom->nlocal + atom->nghost;
            for (int ib = 0; ib < num_b_i; ib++) {
              if (btype_i[ib] <= 1) continue;
              tagint bonded_tag = batom_i[ib];
              int bonded_idx = -1;
              for (int k = 0; k < total_arr; k++) {
                if (tag[k] == bonded_tag) { bonded_idx = k; break; }
              }
              if (bonded_idx < 0) continue;
              // Only intra bonds (same molecule)
              if (atom->molecule[bonded_idx] != mol_i) continue;

              double dx = x[i][0] - x[bonded_idx][0];
              double dy = x[i][1] - x[bonded_idx][1];
              double dz = x[i][2] - x[bonded_idx][2];
              double bd = sqrt(dx*dx + dy*dy + dz*dz);
              if (bd < 1e-6) continue;

              // Soft repulsion: U_rep = epsilon_hb * intra_suppress * exp(-(bd-r_hb)^2/sigma^2)
              // This cancels the intra HB well when inter is formed
              const double eps_rep = 0.08, r_rep = 6.5, sig_rep = 0.8;
              double dr_rep = bd - r_rep;
              double sigma_sq = sig_rep * sig_rep;
              double u_rep = eps_rep * intra_suppress_i * exp(-dr_rep*dr_rep / sigma_sq);
              double f_rep = u_rep * (2.0 * dr_rep / sigma_sq) / bd;

              // Apply to atom i (force away from bonded partner)
              f[i][0] += f_rep * dx;
              f[i][1] += f_rep * dy;
              f[i][2] += f_rep * dz;
              // Apply reaction to bonded_idx if local
              if (bonded_idx < atom->nlocal) {
                f[bonded_idx][0] -= f_rep * dx;
                f[bonded_idx][1] -= f_rep * dy;
                f[bonded_idx][2] -= f_rep * dz;
              }
              // Add to energy (avoid double counting: only when i < bonded_idx)
              if (eflag && i < bonded_idx) {
                evdwl += u_rep;
              }
            }
          }

          // Apply repulsive penalty to intramolecular HB bond partners of j
          if (intra_suppress_j > 1e-6) {
            int num_b_j = atom->num_bond[j];
            int *btype_j = atom->bond_type[j];
            int *batom_j = atom->bond_atom[j];
            tagint *tag = atom->tag;
            int total_arr = atom->nlocal + atom->nghost;
            for (int jb = 0; jb < num_b_j; jb++) {
              if (btype_j[jb] <= 1) continue;
              tagint bonded_tag = batom_j[jb];
              int bonded_idx = -1;
              for (int k = 0; k < total_arr; k++) {
                if (tag[k] == bonded_tag) { bonded_idx = k; break; }
              }
              if (bonded_idx < 0) continue;
              if (atom->molecule[bonded_idx] != mol_j) continue;

              double dx = x[j][0] - x[bonded_idx][0];
              double dy = x[j][1] - x[bonded_idx][1];
              double dz = x[j][2] - x[bonded_idx][2];
              double bd = sqrt(dx*dx + dy*dy + dz*dz);
              if (bd < 1e-6) continue;

              const double eps_rep = 0.08, r_rep = 6.5, sig_rep = 0.8;
              double dr_rep = bd - r_rep;
              double sigma_sq = sig_rep * sig_rep;
              double u_rep = eps_rep * intra_suppress_j * exp(-dr_rep*dr_rep / sigma_sq);
                
                
              double f_rep = u_rep * (2.0 * dr_rep / sigma_sq) / bd;

              f[j][0] += f_rep * dx;
              f[j][1] += f_rep * dy;
              f[j][2] += f_rep * dz;
              if (bonded_idx < atom->nlocal) {
                f[bonded_idx][0] -= f_rep * dx;
                f[bonded_idx][1] -= f_rep * dy;
                f[bonded_idx][2] -= f_rep * dz;
              }
              if (eflag && j < bonded_idx) {
                evdwl += u_rep;
              }
            }
          }

          // Still require different molecules; always enter the block (scale handles suppression)
          if (different_molecules) {
            
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
              
              double r_shifted = r * pow(2.0, 1.0/6.0);  // r_effective = 2^(1/6) * r
              double r_ratio = r_hb[itype][jtype] / r_shifted;
              double barrier = A_barrier[itype][jtype] * pow(r_ratio, 12.0);
              
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
              double dbarrier_dr = -12.0 * barrier * rinv;
              
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
          }  // end different_molecules && !i_has_close_hb_bond && !j_has_close_hb_bond
        }  // end epsilon_hb check

        // ========== APPLY FORCES ==========
        // Radial forces (Morse + Electrostatics): F_radial = fpair * r_hat
        fpair = f_morse + f_elec;
        
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
