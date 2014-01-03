#ifndef _write_PDB_files_h
#define _write_PDB_files_h


#include "md_ld_common_d_dimension_namd.h"
#include "write_pdb_atom_record_v2.h"


void Read_alloc_PDB_info(int *Num_particles, PDB_ATOM **PDB_atoms);

void Write_one_config_to_PDB_file(FR_DAT *fr, PDB_ATOM *PDB_atoms, FILE *out);

void Write_PDB_FILE(FR_DAT *fr, PDB_ATOM *PDB_atoms, char *Filename_PDB);

void Write_config_to_a_PDB_trajectory_file(FR_DAT *fr, PDB_ATOM *PDB_atoms, int config_index, FILE *out);

#endif
