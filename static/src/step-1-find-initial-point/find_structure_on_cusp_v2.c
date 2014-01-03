/* Given two end point strucutres finds a strucutre on the cusp, i.e. a strucutre for
   which energies given by tw osurfaces are same within a tolerance.*/

/* Allows energy offsets for both the end states*/

/* Last updated on May 20, 2013*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "include_all.h"




/* Input file: 'INPUT_FIND_STRUCTURE_ON_CUSP'
 * Format:
 * -----------------------------------------
 *  NUM_PARTICLES (Number of particles in the system %d)
 *  CUT_OFFS (Cutoffs for both the end states. %lf %lf)
 *  FORCE_CONSTANTS (Force constants for both the end states. %lf %lf)
 *  ENERGY_OFFSETS (Energy offsets for both the end states. %lf %lf)
 *  CUSP_TOLERANCE (Energy tolerance for determining strucutre on cusp. %lf)
 *  NUM_IMAGES_INTERPOLATE (Number of images used for interpolation. %d)
 * -----------------------------------------
 * */


int main()
{
  int Num_particles;
  double force_constant_1, force_constant_2;
  double cutoff_1, cutoff_2;
  double energy_offset_1, energy_offset_2;
  double **pair_distances_structure_1, **pair_distances_structure_2;
  double tol;
  int Num_images;
  char Filename[1000], key[1000];
  FR_DAT fr_end_1, fr_end_2, fr_1, fr_2, fr_try_1, fr_try_2, fr_new, fr_cusp;
  FILE *in, *out;
  int i;
  double potential_1, potential_2;

  in = fopen("INPUT_FIND_STRUCTURE_ON_CUSP", "r");
  if(in == NULL)
    {
      printf("File 'INPUT_FIND_STRUCTURE_ON_CUSP' not found.\n");
      exit(1);
    }
  fscanf(in, "%s %d ", key, &Num_particles);
  fscanf(in, "%s %lf %lf ", key, &cutoff_1, &cutoff_2);
  fscanf(in, "%s %lf %lf ", key, &force_constant_1, &force_constant_2);
  fscanf(in, "%s %lf %lf ", key, &energy_offset_1, &energy_offset_2);
  fscanf(in, "%s %lf ", key, &tol);
  fscanf(in, "%s %d ", key, &Num_images); 
  fclose(in);

  pair_distances_structure_1 = (double**)malloc(Num_particles * sizeof(double*));
  pair_distances_structure_2 = (double**)malloc(Num_particles * sizeof(double*));
  for(i = 0 ; i < Num_particles ; i++)
    {
      pair_distances_structure_1[i] = (double*)malloc(Num_particles * sizeof(double));
      pair_distances_structure_2[i] = (double*)malloc(Num_particles * sizeof(double));
    }

  fr_end_1.natoms = Num_particles;
  fr_end_1.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_end_2.natoms = Num_particles ;
  fr_end_2.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_1.natoms = Num_particles;
  fr_1.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_2.natoms = Num_particles;
  fr_2.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_try_1.natoms = Num_particles;
  fr_try_1.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_try_2.natoms = Num_particles;
  fr_try_2.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_new.natoms = Num_particles;
  fr_new.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_new.f = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_cusp.natoms = Num_particles;
  fr_cusp.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_cusp.f = (rvec*)malloc(Num_particles * sizeof(rvec));

  sprintf(Filename, "INPUT_STRUCTURE_1");
  Read_config_simple_3d(Num_particles, &fr_end_1, Filename);
  sprintf(Filename, "INPUT_STRUCTURE_2");
  Read_config_simple_3d(Num_particles, &fr_end_2, Filename);
 
  Calc_pair_distances_both_structures(Num_particles, fr_end_1, fr_end_2, pair_distances_structure_1, pair_distances_structure_2);


  Find_strucutre_on_cusp_offset(Num_particles, &fr_end_1, &fr_end_2, tol, Num_images, pair_distances_structure_1, cutoff_1, force_constant_1, energy_offset_1, pair_distances_structure_2, cutoff_2, force_constant_2, energy_offset_2, &fr_1, &fr_2, &fr_try_1, &fr_try_2, &fr_new, &fr_cusp);

  sprintf(Filename, "initial_struct_on_cusp");
  Write_config_simple(&fr_cusp, Filename);  

  potential_1 = getforces_enm_return_pot_offset(&fr_cusp, pair_distances_structure_1, cutoff_1, force_constant_1, energy_offset_1); 
  potential_2 = getforces_enm_return_pot_offset(&fr_cusp, pair_distances_structure_2, cutoff_2, force_constant_2, energy_offset_2);
  out = fopen("energies_struct_on_cusp", "w");
  fprintf(out, "Energy from surface 1: %.15e\n", potential_1);
  fprintf(out, "Energy from surface 2: %.15e\n", potential_2);
  fclose(out);

  return 0;
}



