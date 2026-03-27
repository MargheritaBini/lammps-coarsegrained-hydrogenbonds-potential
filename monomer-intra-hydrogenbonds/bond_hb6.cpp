// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "bond_hb6.h"
#include <cmath>
#include "atom.h"
#include "neighbor.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondHB62new::BondHB62new(LAMMPS *lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondHB62new::~BondHB62new()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(epsilon);
    memory->destroy(rhb);
    memory->destroy(sigma);
    memory->destroy(A);
  }
}

/* ---------------------------------------------------------------------- */

void BondHB62new::compute(int eflag, int vflag)
{
  int i1, i2, n, type;
  double delx, dely, delz, ebond;
  double rsq, r, rij_hat[3], h1, h2, h0, h, barrier;

  ebond = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    // ========== GEOMETRY CALCULATION ==========
    
    // Compute vector rij = r_i - r_j and its magnitude
    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx * delx + dely * dely + delz * delz;
    if (rsq < 1e-12) continue;
      
    r = sqrt(rsq);
    
    // Unit vector rij_hat = rij / |rij|
    rij_hat[0] = delx / r;
    rij_hat[1] = dely / r;
    rij_hat[2] = delz / r;

    // Check bounds for neighboring atoms (needed for ti and tj)
    if (i1 < 1 || i2 < 1 || i1 >= nlocal - 1 || i2 >= nlocal - 1) continue;
      
    // ========== COMPUTE ti (normal vector at atom i) ==========
    
    double vi1[3], vi2[3], ti[3], ti_unit[3];
    double ti_norm;

    // vi1 = r_i - r_{i-1}
    // vi2 = r_{i+1} - r_i
    for (int d = 0; d < 3; d++) {
      vi1[d] = x[i1][d] - x[i1 - 1][d];
      vi2[d] = x[i1 + 1][d] - x[i1][d];
    }

    // Cross product ti = vi1 × vi2
    ti[0] = vi1[1] * vi2[2] - vi1[2] * vi2[1];
    ti[1] = vi1[2] * vi2[0] - vi1[0] * vi2[2];
    ti[2] = vi1[0] * vi2[1] - vi1[1] * vi2[0];

    // Normalize ti
    ti_norm = sqrt(ti[0] * ti[0] + ti[1] * ti[1] + ti[2] * ti[2]);

    if (ti_norm < 1e-12) {
      // Handle degenerate case - set to zero
      ti_unit[0] = ti_unit[1] = ti_unit[2] = 0.0;
    } else {
      for (int d = 0; d < 3; d++) {
        ti_unit[d] = ti[d] / ti_norm;
      }
    }

    // ========== COMPUTE tj (normal vector at atom j) ==========
    
    double vj1[3], vj2[3], tj[3], tj_unit[3];
    double tj_norm;

    // vj1 = r_j - r_{j-1}
    // vj2 = r_{j+1} - r_j
    for (int d = 0; d < 3; d++) {
      vj1[d] = x[i2][d] - x[i2 - 1][d];
      vj2[d] = x[i2 + 1][d] - x[i2][d];
    }

    // Cross product tj = vj1 × vj2
    tj[0] = vj1[1] * vj2[2] - vj1[2] * vj2[1];
    tj[1] = vj1[2] * vj2[0] - vj1[0] * vj2[2];
    tj[2] = vj1[0] * vj2[1] - vj1[1] * vj2[0];

    // Normalize tj
    tj_norm = sqrt(tj[0] * tj[0] + tj[1] * tj[1] + tj[2] * tj[2]);

    if (tj_norm < 1e-12) {
      // Handle degenerate case - set to zero
      tj_unit[0] = tj_unit[1] = tj_unit[2] = 0.0;
    } else {
      for (int d = 0; d < 3; d++) {
        tj_unit[d] = tj[d] / tj_norm;
      }
    }

    // ========== COMPUTE ANGULAR TERMS ==========
    
    // alpha_i = rij_hat · ti_unit
    double alpha_i = 0.0;
    for (int d = 0; d < 3; d++) {
      alpha_i += rij_hat[d] * ti_unit[d];
    }

    // alpha_j = rij_hat · tj_unit
    double alpha_j = 0.0;
    for (int d = 0; d < 3; d++) {
      alpha_j += rij_hat[d] * tj_unit[d];
    }

    // ========== COMPUTE ENERGY TERMS ==========
    
    double sigma_sq = sigma[type] * sigma[type];
    
    // H0 = -epsilon * exp(-(r - r_hb)^2 / sigma^2)
    h0 = -epsilon[type] * exp(-pow(r - rhb[type], 2) / sigma_sq);
    
    // H1 = exp((|alpha_i| - 1) / sigma^2)
    h1 = exp((fabs(alpha_i) - 1.0) / sigma_sq);
    
    // H2 = exp((|alpha_j| - 1) / sigma^2)
    h2 = exp((fabs(alpha_j) - 1.0) / sigma_sq);
    
    // Barrier term = A * (r_hb / r)^12
    double r_eff = rhb[type] / pow(2.0, 1.0/6.0);  // r_hb / 2^(1/6)
    barrier = A[type] * pow(r_eff / r, 12.0);
    
    // Total energy U = H0 * H1 * H2 + barrier
    h = h0 * h1 * h2 + barrier;

    // ========== FORCE CALCULATION ==========
    // Based on frozen backbone approximation where ti_unit and tj_unit are constant
    
    // --- Step 1: Derivative of H0 with respect to r ---
    // dH0/dr = H0 * (-2(r - r_hb) / sigma^2)
    double dH0_dr = h0 * (-2.0) * (r - rhb[type]) / sigma_sq;
    
    // --- Step 2: Derivative of barrier with respect to r ---
    // dU_barrier/dr = -12 * barrier / r
    double dbarrier_dr = -12.0 * barrier / r;
    
    // --- Step 3: Sign function for alpha_i and alpha_j ---
    double sgn_alpha_i, sgn_alpha_j;
    
    if (fabs(alpha_i) < 1e-12) {
      sgn_alpha_i = 0.0;
    } else {
      sgn_alpha_i = (alpha_i > 0.0) ? 1.0 : -1.0;
    }
    
    if (fabs(alpha_j) < 1e-12) {
      sgn_alpha_j = 0.0;
    } else {
      sgn_alpha_j = (alpha_j > 0.0) ? 1.0 : -1.0;
    }
    
    // --- Step 4: Compute gradient vector components for H1 ---
    // dH1/dr_i = (H1 * sgn(alpha_i) / (sigma^2 * r)) * (ti_unit - alpha_i * rij_hat)
    double dH1_coeff = (h1 * sgn_alpha_i) / (sigma_sq * r);
    double dH1_dri[3];
    for (int d = 0; d < 3; d++) {
      dH1_dri[d] = dH1_coeff * (ti_unit[d] - alpha_i * rij_hat[d]);
    }
    
    // --- Step 5: Compute gradient vector components for H2 ---
    // dH2/dr_i = (H2 * sgn(alpha_j) / (sigma^2 * r)) * (tj_unit - alpha_j * rij_hat)
    double dH2_coeff = (h2 * sgn_alpha_j) / (sigma_sq * r);
    double dH2_dri[3];
    for (int d = 0; d < 3; d++) {
      dH2_dri[d] = dH2_coeff * (tj_unit[d] - alpha_j * rij_hat[d]);
    }
    
    // --- Step 6: Total gradient dU/dr_i using product rule ---
    // dU/dr_i = (dH0/dr * rij_hat) * H1 * H2
    //         + H0 * (dH1/dr_i) * H2
    //         + H0 * H1 * (dH2/dr_i)
    //         + (dbarrier/dr * rij_hat)
    
    double dU_dri[3];
    for (int d = 0; d < 3; d++) {
      dU_dri[d] = (dH0_dr * rij_hat[d]) * h1 * h2           // H0 term
                  + h0 * dH1_dri[d] * h2                     // H1 term
                  + h0 * h1 * dH2_dri[d]                     // H2 term
                  + dbarrier_dr * rij_hat[d];                // barrier term
    }
    
    // --- Step 7: Compute forces ---
    // Force on atom i1: F_i = -dU/dr_i
    // Force on atom i2: F_j = +dU/dr_i (Newton's third law)
    
    double fi[3], fj[3];
    for (int d = 0; d < 3; d++) {
      fi[d] = -dU_dri[d];
      fj[d] = +dU_dri[d];
    }

    // ========== APPLY FORCES ==========
    
    if (newton_bond || i1 < nlocal) {
      f[i1][0] += fi[0];
      f[i1][1] += fi[1];
      f[i1][2] += fi[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += fj[0];
      f[i2][1] += fj[1];
      f[i2][2] += fj[2];
    }

    // ========== ENERGY AND VIRIAL ==========
    
    if (eflag) ebond = h;

    if (evflag) {
      // For virial calculation, we need the scalar force magnitude
      // along the bond direction: fbond = -F_i · rij_hat
      double fbond = -(fi[0] * delx + fi[1] * dely + fi[2] * delz) / r;
      ev_tally(i1, i2, nlocal, newton_bond, ebond, fbond, delx, dely, delz);
    }
  }
}

/* ---------------------------------------------------------------------- */

void BondHB62new::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(epsilon, n + 1, "bond:epsilon");
  memory->create(rhb, n + 1, "bond:rhb");
  memory->create(sigma, n + 1, "bond:sigma");
  memory->create(A, n + 1, "bond:A");
  memory->create(setflag, n + 1, "bond:setflag");

  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ---------------------------------------------------------------------- */

void BondHB62new::coeff(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR, "Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double epsilon_one = utils::numeric(FLERR, arg[1], false, lmp);
  double rhb_one = utils::numeric(FLERR, arg[2], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[3], false, lmp);
  double A_one = utils::numeric(FLERR, arg[4], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
        epsilon[i] = epsilon_one;
        rhb[i] = rhb_one;
        sigma[i] = sigma_one;
        A[i] = A_one;
        setflag[i] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
}


double BondHB62new::equilibrium_distance(int i)
{
  return rhb[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondHB62new::write_restart(FILE *fp)
{
  fwrite(&epsilon[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&rhb[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&sigma[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&A[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondHB62new::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &epsilon[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &rhb[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &sigma[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &A[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&epsilon[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&rhb[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&sigma[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&A[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondHB62new::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++) {
    fprintf(fp, "%d %g %g %g %g\n", i, epsilon[i], rhb[i], sigma[i], A[i]);
  }
}

/* ---------------------------------------------------------------------- */

double BondHB62new::single(int type, double rsq, int i, int j, double &fforce)
{
  // Declare necessary variables
  double delx, dely, delz, rij_hat[3], r, h0, h1, h2, h, barrier;
  
  // Get positions of atoms i and j
  double **x = atom->x;

  // Compute vector rij and its magnitude
  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];

  // Calculate distance
  rsq = delx * delx + dely * dely + delz * delz;
  r = sqrt(rsq);
  
  // Unit vector
  rij_hat[0] = delx / r;
  rij_hat[1] = dely / r;
  rij_hat[2] = delz / r;

  // Compute the exponential decay term (H0)
  double sigma_sq = sigma[type] * sigma[type];
  h0 = -epsilon[type] * exp(-pow(r - rhb[type], 2) / sigma_sq);

  // Compute angular dependence (H1 and H2)
  double vi1[3], vi2[3], ti[3], ti_unit[3];
  double vj1[3], vj2[3], tj[3], tj_unit[3];
  double ti_norm, tj_norm;

  // Calculate cross products for angular terms
  for (int d = 0; d < 3; d++) {
    vi1[d] = x[i][d] - x[i - 1][d];
    vi2[d] = x[i + 1][d] - x[i][d];
  }

  ti[0] = vi1[1] * vi2[2] - vi1[2] * vi2[1];
  ti[1] = vi1[2] * vi2[0] - vi1[0] * vi2[2];
  ti[2] = vi1[0] * vi2[1] - vi1[1] * vi2[0];

  // Normalize the vector
  ti_norm = sqrt(ti[0] * ti[0] + ti[1] * ti[1] + ti[2] * ti[2]);
  
  if (ti_norm < 1e-12) {
    ti_unit[0] = ti_unit[1] = ti_unit[2] = 0.0;
  } else {
    for (int d = 0; d < 3; d++) ti_unit[d] = ti[d] / ti_norm;
  }

  // Repeat similar calculations for the second atom (vj and tj)
  for (int d = 0; d < 3; d++) {
    vj1[d] = x[j][d] - x[j - 1][d];
    vj2[d] = x[j + 1][d] - x[j][d];
  }

  tj[0] = vj1[1] * vj2[2] - vj1[2] * vj2[1];
  tj[1] = vj1[2] * vj2[0] - vj1[0] * vj2[2];
  tj[2] = vj1[0] * vj2[1] - vj1[1] * vj2[0];

  tj_norm = sqrt(tj[0] * tj[0] + tj[1] * tj[1] + tj[2] * tj[2]);
  
  if (tj_norm < 1e-12) {
    tj_unit[0] = tj_unit[1] = tj_unit[2] = 0.0;
  } else {
    for (int d = 0; d < 3; d++) tj_unit[d] = tj[d] / tj_norm;
  }

  // Calculate dot products for angular terms
  double alpha_i = 0.0, alpha_j = 0.0;
  for (int d = 0; d < 3; d++) {
    alpha_i += rij_hat[d] * ti_unit[d];
    alpha_j += rij_hat[d] * tj_unit[d];
  }

  // Compute H1 and H2
  h1 = exp((fabs(alpha_i) - 1.0) / sigma_sq);
  h2 = exp((fabs(alpha_j) - 1.0) / sigma_sq);

  // Compute barrier
  double r_eff = rhb[type] / pow(2.0, 1.0/6.0);
  barrier = A[type] * pow(r_eff / r, 12.0);
  
  // Total energy
  h = h0 * h1 * h2 + barrier;
  
  // ========== FORCE CALCULATION ==========
  
  // Derivative of H0
  double dH0_dr = h0 * (-2.0) * (r - rhb[type]) / sigma_sq;
  
  // Derivative of barrier
  double dbarrier_dr = -12.0 * barrier / r;
  
  // Sign functions
  double sgn_alpha_i = (fabs(alpha_i) < 1e-12) ? 0.0 : ((alpha_i > 0.0) ? 1.0 : -1.0);
  double sgn_alpha_j = (fabs(alpha_j) < 1e-12) ? 0.0 : ((alpha_j > 0.0) ? 1.0 : -1.0);
  
  // Gradient components for H1
  double dH1_coeff = (h1 * sgn_alpha_i) / (sigma_sq * r);
  double dH1_dri[3];
  for (int d = 0; d < 3; d++) {
    dH1_dri[d] = dH1_coeff * (ti_unit[d] - alpha_i * rij_hat[d]);
  }
  
  // Gradient components for H2
  double dH2_coeff = (h2 * sgn_alpha_j) / (sigma_sq * r);
  double dH2_dri[3];
  for (int d = 0; d < 3; d++) {
    dH2_dri[d] = dH2_coeff * (tj_unit[d] - alpha_j * rij_hat[d]);
  }
  
  // Total gradient
  double dU_dri[3];
  for (int d = 0; d < 3; d++) {
    dU_dri[d] = (dH0_dr * rij_hat[d]) * h1 * h2
                + h0 * dH1_dri[d] * h2
                + h0 * h1 * dH2_dri[d]
                + dbarrier_dr * rij_hat[d];
  }
  
  // Force magnitude along bond direction
  // fbond = -F_i · rij_hat = dU_dri · rij_hat
  double fbond = 0.0;
  for (int d = 0; d < 3; d++) {
    fbond += dU_dri[d] * rij_hat[d];
  }
  
  fforce = fbond;

  return h;  // Return the computed bond energy
}
