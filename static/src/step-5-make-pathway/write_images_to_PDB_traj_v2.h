#ifndef _write_images_to_PDB_traj_v2_h
#define _write_images_to_PDB_traj_v2_h

#include "md_ld_common_d_dimension_namd.h"
#include "read_configuration.h"
#include "write_pdb_atom_record_v2.h"
#include "write_pdb_atom_record_v2.h"
#include "write_PDB_files.h"

void Make_PDB_traj_from_images_v2(int Num_images, int Num_atoms, char *Filename_PDB_traj);

#endif
