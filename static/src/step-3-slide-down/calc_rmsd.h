#ifndef _calc_rmsd_h
#define _calc_rmsd_h


#include "md_ld_common_d_dimension_namd.h"
#include "read_configuration.h"

double Calc_rmsd_one_structure(FR_DAT *fr, FR_DAT *fr_fixed_state);

void Calc_rmsd_all_structure(int Num_atoms, int Num_images, char *Filename_base, char *Filename_fixed_base, FILE *out_rmsd_structure, long int iteration_number);

double Calc_rmsd_one_image_of_string(int Image_dimension, double *Image_position, double *Image_position_fixed);

void Calc_rmsd_string(int Num_images, int Image_dimension, char *Filename_base, char *Filename_fixed_base, FILE *out_all_rmsd, FILE *out_average_rmsd, long int iteration_number);

#endif
