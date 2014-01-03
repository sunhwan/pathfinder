#ifndef _force_pot_simple_enm_h
#define _force_pot_simple_enm_h

#include "md_ld_common_d_dimension_namd.h"

void Calc_simple_distance(int dimension, double *point_1, double *point_2, double *distance, double *distance_squared, double *distance_vector);

void Calc_pair_distances_both_structures(int Num_particles, FR_DAT fr_structure_1, FR_DAT fr_structure_2, double **pair_distances_structure_1, double **pair_distances_structure_2);

void getforces_enm(FR_DAT *fr, double **pair_distances_structure, double cutoff, double force_constant);

double getforces_enm_return_pot(FR_DAT *fr, double **pair_distances_structure, double cutoff, double force_constant);

void getforces_enm_offset(FR_DAT *fr, double **pair_distances_structure, double cutoff, double force_constant, double energy_offset);

double getforces_enm_return_pot_offset(FR_DAT *fr, double **pair_distances_structure, double cutoff, double force_constant, double energy_offset);

#endif
