#ifndef _wrtie_configuration_h
#define _wrtie_configuration_h

#include "md_ld_common_d_dimension_namd.h"

void Write_config(FR_DAT *fr, char *output_filename);

void Write_config_to_trajectory_file(FR_DAT *fr, FILE *output_file_pointer);

void Write_config_to_trajectory_file_only_coordinates(FR_DAT *fr, FILE *output_file_pointer);

void Write_config_simple(FR_DAT *fr, char *output_filename);

void Write_config_vmd_xyz(FR_DAT *fr, char *output_filename);

#endif
