#ifndef _calc_energy_pathway_v2_h
#define _calc_energy_pathway_v2_h

#include "calc_rmsd.h"
#include "force_pot_simple_enm.h"
#include "md_ld_common_d_dimension_namd.h"
#include "read_configuration.h"


void Calc_energies_along_pathway_offset(int Num_particles, FR_DAT *fr_ref_1, FR_DAT *fr_ref_2, double force_constant_1, double force_constant_2, double cutoff_1, double cutoff_2, double energy_offset_1, double energy_offset_2, double **pair_distances_structure_1, double **pair_distances_structure_2, int Num_images, int index_of_ts_struct);

#endif
