/* Calculate energy along the pathway.*/

/* Last modified on May 20, 2013*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "calc_energy_pathway_v2.h"




/* Energies of all the images along a pathway. 
   'index_of_ts_struct' is the index of the transition state structure, index is 1-based.*/
void Calc_energies_along_pathway_offset(int Num_particles, FR_DAT *fr_ref_1, FR_DAT *fr_ref_2, double force_constant_1, double force_constant_2, double cutoff_1, double cutoff_2, double energy_offset_1, double energy_offset_2, double **pair_distances_structure_1, double **pair_distances_structure_2, int Num_images, int index_of_ts_struct)
{
  FR_DAT fr;
  char Filename[1000];
  int alpha;
  FILE *out;
  double energy_1, energy_2;

  fr.natoms = Num_particles;
  fr.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr.f = (rvec*)malloc(Num_particles * sizeof(rvec));

  out = fopen("energies_pathway", "w");
  for(alpha = 0 ; alpha < Num_images ; alpha++)
    {
      if((alpha+1) < index_of_ts_struct) /* Surface 1*/
	{
	  sprintf(Filename, "COORDS_IMAGE_%d", (alpha+1));
	  Read_config_simple_3d(Num_particles, &fr, Filename);
	  getforces_enm_offset(&fr, pair_distances_structure_1, cutoff_1, force_constant_1, energy_offset_1);
	  energy_1 = fr.U;
	  fprintf(out, "%d  %.15e\n", (alpha+1), energy_1);
	}
      else if((alpha+1) == index_of_ts_struct) /* TS*/
	{
	  sprintf(Filename, "COORDS_IMAGE_%d", (alpha+1));
	  Read_config_simple_3d(Num_particles, &fr, Filename);
	  /* Calcualte using surface 1*/
	  getforces_enm_offset(&fr, pair_distances_structure_1, cutoff_1, force_constant_1, energy_offset_1);
	  energy_1 = fr.U;
	  /* Calcualte using surface 2*/
	  getforces_enm_offset(&fr, pair_distances_structure_2, cutoff_2, force_constant_2, energy_offset_2);
	  energy_2 = fr.U;
	  fprintf(out, "%d  %.15e %.15e\n", (alpha+1), energy_1, energy_2);
	}
      else
	{
	  sprintf(Filename, "COORDS_IMAGE_%d", (alpha+1));
	  Read_config_simple_3d(Num_particles, &fr, Filename);
	  getforces_enm_offset(&fr, pair_distances_structure_2, cutoff_2, force_constant_2, energy_offset_2);
	  energy_2 = fr.U;
	  fprintf(out, "%d  %.15e\n", (alpha+1), energy_2);
	}
    }
  fclose(out);
}







