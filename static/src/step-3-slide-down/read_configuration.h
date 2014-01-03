#ifndef _read_configuration_h
#define _read_configuration_h

#include "md_ld_common_d_dimension_namd.h"

void Read_config(FR_DAT *fr, char *input_filename, int Num_particles_from_top);

void Read_config_alloc(FR_DAT *fr, char *input_filename);

void Read_alloc_config_simple_3d(int Num_atoms, FR_DAT *fr, char *input_filename);

void Read_config_simple_3d(int Num_atoms, FR_DAT *fr, char *input_filename);

void Read_config_vmd_xyz(FR_DAT *fr, char *input_filename);

void Read_alloc_config_vmd_xyz(FR_DAT *fr, char *input_filename);

#endif
