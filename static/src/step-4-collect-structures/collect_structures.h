#ifndef _collect_structures_h
#define _collect_structures_h

#include "md_ld_common_d_dimension_namd.h"
#include "read_configuration.h"
#include "wrtie_configuration.h"

void Collect_structures(int Num_particles, int Num_structures_1, int Num_structures_2, int *Num_struct_in_pathway);

void Collect_structures_one_ts(int Num_particles, FR_DAT *fr_ts, int Num_structures_1, int Num_structures_2, int *Num_struct_in_pathway, int *index_of_ts_struct);

#endif
