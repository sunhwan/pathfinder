#ifndef _slide_down_v3_h
#define _slide_down_v3_h

#include "align_two_structures_Lei.h"
#include "calc_rmsd.h"
#include "force_pot_simple_enm.h"
#include "md_ld_common_d_dimension_namd.h"
#include "wrtie_configuration.h"



int Check_step_size_offset(int Num_particles, double force_constant, double cutoff, double energy_offset, double **pair_distances_structure, double *Reference_x, double *Reference_y, double *Reference_z, double rmsd_pathway, double delta_rmsd_zero, FR_DAT *fr_start, double step_size);

void Slide_down_single_step_size_offset(int Num_particles, double force_constant, double cutoff, double energy_offset, double **pair_distances_structure, long int Max_num_iterations, double *Reference_x, double *Reference_y, double *Reference_z, int struct_index, double rmsd_pathway, double delta_rmsd_zero, FR_DAT *fr_ref, FR_DAT *fr_start, double step_size, double energy_ref, double energy_from_ref_tol, int *Num_struct_written, FR_DAT *fr_end, long int *iteration_number_end, double step_size_reduction_factor);

#endif
