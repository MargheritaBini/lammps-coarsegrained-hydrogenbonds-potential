/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   Simplified HB/Elec/Morse pair style
   
   Simplified parameter format:
   pair_coeff i j D0 alpha r0 lambda epsilon_hb r_hb sigma A
   
   where:
   - lambda = unified charge parameter (qi * qj product)
   - sigma = unified width (same for distance and angular)
   - No separate donor/acceptor flags (symmetric potential)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(hb/nonbonded,PairHBNonbonded);
// clang-format on
#else

#ifndef LMP_PAIR_HB_NONBONDED_H
#define LMP_PAIR_HB_NONBONDED_H

#include "pair.h"

namespace LAMMPS_NS {

class PairHBNonbonded : public Pair {
 public:
  PairHBNonbonded(class LAMMPS *);
  virtual ~PairHBNonbonded();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);

 protected:
  // Morse potential parameters
  double **d0, **alpha, **r0;
  
  // Unified electrostatic parameter (qi*qj product)
  double **lambda;
  
  // Hydrogen bond parameters (simplified)
  double **epsilon_hb;
  double **r_hb;
  double **sigma;        // Unified width parameter
  double **A_barrier;
  
  // Global parameters
  double kappa;          // Screening parameter
  double dielectric;     // Dielectric constant
  double qqrd2e;         // Charge conversion factor
  
  // Control flags
  int use_hb_exclusion;  // Flag to exclude electrostatics when HB active
  
  double cut_global;
  double **cut;
  double **offset;

  void allocate();
  void compute_normal_vector(double **, int, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
