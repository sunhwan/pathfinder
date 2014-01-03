/* With the TS structures, end strucutres and the strucutres obtained
   by two separate steepest descent minimization this program performs the 
   following taks
   1. Collect structures in the following order to make the pathway
      End point 1 -> Strucutres on the surface 1 in reverse order of their generation
      -> TS -> Structures on surface 2 in the order they were generated -> End point 2.
   2. Calculate energy along the pathway.
*/

/* Allows energy offsets for both the end states*/

/* Last modified on May 20, 2013*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "include_all.h"


/* Input filename: 'INPUT_COLLECT_ENERGY'
 * Format:
 * ---------------------------
 *  NUM_PARTICLES (Number of particles in the system. %d)
 *  CUT_OFFS (Cutoffs for both the end states. %lf %lf)
 *  FORCE_CONSTANTS (Force constants for both the end states. %lf %lf)
 *  ENERGY_OFFSETS (Energy offsets for both the end states. %lf %lf)
 * ---------------------------
 * */


void Read_inputs(int *Num_particles, double *cutoff_1, double *cutoff_2, double *force_constant_1, double *force_constant_2, double *energy_offset_1, double *energy_offset_2)
{
  FILE *in;
  char key[1000];

  in = fopen("INPUT_COLLECT_ENERGY", "r");
  if(in == NULL)
    {
      printf("File 'INPUT_COLLECT_ENERGY' not found\n");
      exit(1);
    }
  fscanf(in, "%s  %d ", key, Num_particles);
  fscanf(in, "%s %lf %lf ", key, cutoff_1, cutoff_2);
  fscanf(in, "%s %lf %lf ", key, force_constant_1, force_constant_2);
  fscanf(in, "%s %lf %lf ", key, energy_offset_1, energy_offset_2);
  fclose(in);
}


int main()
{
  int Num_particles;
  double cutoff_1, cutoff_2;
  double force_constant_1, force_constant_2;
  double energy_offset_1, energy_offset_2;
  FR_DAT fr_ref_1, fr_ref_2, fr_ts, fr_1, fr_2;
  double **pair_distances_structure_1, **pair_distances_structure_2;
  char Filename[1000];
  int i;
  FILE *in;
  int Num_struct_written_surface_1, Num_struct_written_surface_2;
  int Num_struct_in_pathway;
  int index_of_ts_struct;
  double rmsd_consecutive;
  FILE *out_rmsd;

  /****************************************************************************************/

  Read_inputs(&Num_particles, &cutoff_1, &cutoff_2, &force_constant_1, &force_constant_2, &energy_offset_1, &energy_offset_2);

  in = fopen("num_structures_written_1", "r");
  if(in == NULL)
    {
      printf("File 'num_structures_written_1' not found\n");
      exit(1);
    }
  fscanf(in, "%d ", &Num_struct_written_surface_1);
  fclose(in);
  in = fopen("num_structures_written_2", "r");
  if(in == NULL)
    {
      printf("File 'num_structures_written_2' not found\n");
      exit(1);
    }
  fscanf(in, "%d ", &Num_struct_written_surface_2);
  fclose(in);

  pair_distances_structure_1 = (double**)malloc(Num_particles * sizeof(double*));
  pair_distances_structure_2 = (double**)malloc(Num_particles * sizeof(double*));
  for(i = 0 ; i < Num_particles ; i++)
    {
      pair_distances_structure_1[i] = (double*)malloc(Num_particles * sizeof(double));
      pair_distances_structure_2[i] = (double*)malloc(Num_particles * sizeof(double));
    }

  fr_ref_1.natoms = Num_particles;
  fr_ref_1.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_ref_2.natoms = Num_particles;
  fr_ref_2.x = (rvec*)malloc(Num_particles * sizeof(rvec));

  fr_ts.natoms = Num_particles;
  fr_ts.x = (rvec*)malloc(Num_particles * sizeof(rvec));

  fr_1.natoms = Num_particles;
  fr_1.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_2.natoms = Num_particles;
  fr_2.x = (rvec*)malloc(Num_particles * sizeof(rvec));


  sprintf(Filename, "INPUT_STRUCTURE_1");
  Read_config_simple_3d(Num_particles, &fr_ref_1, Filename);
  sprintf(Filename, "INPUT_STRUCTURE_2");
  Read_config_simple_3d(Num_particles, &fr_ref_2, Filename);
 
  Calc_pair_distances_both_structures(Num_particles, fr_ref_1, fr_ref_2, pair_distances_structure_1, pair_distances_structure_2);

  sprintf(Filename, "minimized_struct_on_cusp");
  Read_config_simple_3d(Num_particles, &fr_ts, Filename);


  /* Collect strucutres*/
  Collect_structures_one_ts(Num_particles, &fr_ts, Num_struct_written_surface_1, Num_struct_written_surface_2, &Num_struct_in_pathway, &index_of_ts_struct);


  /* Calculate energy*/
  Calc_energies_along_pathway_offset(Num_particles, &fr_ref_1, &fr_ref_2, force_constant_1, force_constant_2, cutoff_1, cutoff_2, energy_offset_1, energy_offset_2, pair_distances_structure_1, pair_distances_structure_2, Num_struct_in_pathway, index_of_ts_struct);

  /* RMSD of consecutive strucutres along the pathway*/
  out_rmsd = fopen("consecutive_rmsd", "w");
  for(i = 0 ; i < (Num_struct_in_pathway-1) ; i++)
    {
      sprintf(Filename, "COORDS_IMAGE_%d", (i+1));
      Read_config_simple_3d(Num_particles, &fr_1, Filename);
      sprintf(Filename, "COORDS_IMAGE_%d", (i+2));
      Read_config_simple_3d(Num_particles, &fr_2, Filename);
      rmsd_consecutive = Calc_rmsd_one_structure(&fr_1, &fr_2); 
      fprintf(out_rmsd, "%.8f\n", rmsd_consecutive);
    }
  fclose(out_rmsd);

  return 0;
}












