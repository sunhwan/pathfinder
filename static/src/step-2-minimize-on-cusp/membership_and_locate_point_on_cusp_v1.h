#ifndef _membership_and_locate_point_on_cusp_v1_h
#define _membership_and_locate_point_on_cusp_v1_h


#include "md_ld_common_d_dimension_namd.h"
#include "force_pot_simple_enm.h"


int Find_membership_offset(FR_DAT *fr, double **pair_distances_structure_1, double cutoff_1, double force_constant_1, double energy_offset_1, double **pair_distances_structure_2, double cutoff_2, double force_constant_2, double energy_offset_2, double tol);

void Gen_struct_from_param(int Num_atoms, FR_DAT *fr_1, FR_DAT *fr_2, double t, FR_DAT *fr_new);

void Copy_config(int Num_atoms, FR_DAT *fr_old, FR_DAT *fr_new);

void Find_strucutre_on_cusp_offset(int Num_atoms, FR_DAT *fr_end_1, FR_DAT *fr_end_2, double tol, int Num_images, double **pair_distances_structure_1, double cutoff_1, double force_constant_1, double energy_offset_1, double **pair_distances_structure_2, double cutoff_2, double force_constant_2, double energy_offset_2, FR_DAT *fr_1, FR_DAT *fr_2, FR_DAT *fr_try_1, FR_DAT *fr_try_2, FR_DAT *fr_new, FR_DAT *fr_cusp);

#endif
