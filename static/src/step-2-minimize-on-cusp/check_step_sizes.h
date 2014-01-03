#ifndef _check_step_sizes_h
#define _check_step_sizes_h

#include "align_two_structures_Lei.h"
#include "calc_rmsd.h"
#include "force_pot_simple_enm.h"
#include "md_ld_common_d_dimension_namd.h"
#include "membership_and_locate_point_on_cusp_v1.h"
#include "read_configuration.h"


int Is_energy_decreasing(int Num_points, double *energy_values);


void Choose_step_size_offset(int Num_particles, double force_constant_1, double force_constant_2, double cutoff_1, double cutoff_2, double energy_offset_1, double energy_offset_2, double **pair_distances_structure_1, double **pair_distances_structure_2, double tol, int Num_images, FR_DAT *fr_end_1, FR_DAT *fr_end_2, FR_DAT *fr_1, FR_DAT *fr_2, FR_DAT *fr_try_1, FR_DAT *fr_try_2, FR_DAT *fr_new, FR_DAT *fr_cusp, double *Reference_x, double *Reference_y, double *Reference_z, int Num_points, double step_size_reduction_factor_1, double step_size_reduction_factor_2, double *step_size_1, double *step_size_2);
#endif
