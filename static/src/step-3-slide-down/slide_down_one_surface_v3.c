/* Starting from the energy minimized strucutre on the cusp region (the transition
   state) this program performs the following task
   1. Slides down on one well.
 */

/* Chooses step-size adaptively if the input value does not work.*/

/* Allows energy offsets for both the end states*/

/* Last modified on May 20, 2013*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "include_all.h"


/* Input filename: 'INPUT_SLIDE_ONE_SURFACE'
 * Format:
 * ---------------------------
 *  NUM_PARTICLES (Number of particles in the system. %d)
 *  SURFACE_INDEX (1 or 2. %d)
 *  CUT_OFF (Cutoff. %lf)
 *  FORCE_CONSTANT (Force constant. %lf)
 *  ENERGY_OFFSET (Energy offset for the end state. %lf)
 *  STEP_SIZE (Step size for steepest descent minimization. %lf)
 *  STEP_SIZE_REDUCTION_FACTOR (Factor which step-size is decreased. %lf)
 *  RMSD_PATHWAY (RMSD of two successive strucutres in the pathway. %lf)
 *  RMSD_DIFF_TOL (Two rmsd values are same if they differ by less than this. %lf)
 *  ENERGY_FROM_REFERENCE_TOL (Energy from the reference end strucutre tolerance. %lf)
 *  MAX_NUM_ITER (Maximum number of iterations. %lu)
 * ---------------------------
 * */

void Read_input(int *Num_particles, int *struct_index, double *cutoff, double *force_constant, double *energy_offset, double *step_size, double *step_size_reduction_factor, double *rmsd_pathway, double *delta_rmsd_zero, double *energy_from_ref_tol, long int *Max_num_iterations)
{
  FILE *in;
  char key[1000];
  in = fopen("INPUT_SLIDE_ONE_SURFACE", "r");
  if(in == NULL)
    {
      printf("File 'INPUT_SLIDE_ONE_SURFACE' not found\n");
      exit(1);
    }
  fscanf(in, "%s %d", key, Num_particles);
  fscanf(in, "%s %d", key, struct_index);
  fscanf(in, "%s %lf", key, cutoff);
  fscanf(in, "%s %lf", key, force_constant);
  fscanf(in, "%s %lf", key, energy_offset);
  fscanf(in, "%s %lf", key, step_size);
  fscanf(in, "%s %lf", key, step_size_reduction_factor);
  fscanf(in, "%s %lf", key, rmsd_pathway);
  fscanf(in, "%s %lf", key, delta_rmsd_zero);
  fscanf(in, "%s %lf", key, energy_from_ref_tol);
  fscanf(in, "%s %lu", key, Max_num_iterations);
  fclose(in);
}


/* Calculate pair distances for one end strucutre*/
void Calc_pair_distances_one_structure(int Num_particles, FR_DAT fr_structure, double **pair_distances_structure)
{
  int i, j;
  int dimension = 3;
  double distance;
  double distance_squared;
  double distance_vector[3];

  /* Set diagonal elements to -1.0*/
  for(i = 0 ; i < Num_particles ; i++)
    {
      pair_distances_structure[i][i] = -1.0;
    }

  for(i = 0 ; i < (Num_particles - 1) ; i++)
    {
      for(j = (i + 1) ; j < Num_particles ; j++)
        {
          /* Structure*/
          Calc_simple_distance(dimension, fr_structure.x[i], fr_structure.x[j], &distance, &distance_squared, distance_vector);
          pair_distances_structure[i][j] = distance;
          pair_distances_structure[j][i] = distance;
        }
    }
}




int main()
{
  int Num_particles;
  int struct_index;
  double cutoff;
  double force_constant;
  double energy_offset;
  double step_size;
  double step_size_reduction_factor;
  double rmsd_pathway;
  double delta_rmsd_zero;
  double energy_from_ref_tol;
  long int Max_num_iterations;

  FR_DAT fr_ref, fr_ref_for_alignment, fr_ts, fr_start, fr_end; 
  double **pair_distances_structure;
  char Filename[1000];
  int i;
  double *Reference_x, *Reference_y, *Reference_z;
  int Num_struct_written;
  FILE *out_num, *out;
  double energy_ref, rmsd_from_ref_final;
  long int iteration_number_end;
  int is_step_size_adequate;

  /********************************************************************************/

  Read_input(&Num_particles, &struct_index, &cutoff, &force_constant, &energy_offset, &step_size, &step_size_reduction_factor, &rmsd_pathway, &delta_rmsd_zero, &energy_from_ref_tol, &Max_num_iterations);

  pair_distances_structure = (double**)malloc(Num_particles * sizeof(double*));
  for(i = 0 ; i < Num_particles ; i++)
    {
      pair_distances_structure[i] = (double*)malloc(Num_particles * sizeof(double));
    }

  fr_ref.natoms = Num_particles;
  fr_ref.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_ref.f = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_ref_for_alignment.natoms = Num_particles;
  fr_ref_for_alignment.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_ts.natoms = Num_particles;
  fr_ts.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_start.natoms = Num_particles;
  fr_start.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_end.natoms = Num_particles;
  fr_end.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_end.f = (rvec*)malloc(Num_particles * sizeof(rvec));

  /*Lei's code uses array whose indices start from 1.*/
  Reference_x = (double*)malloc((Num_particles + 1) * sizeof(double));
  Reference_y = (double*)malloc((Num_particles + 1) * sizeof(double));
  Reference_z = (double*)malloc((Num_particles + 1) * sizeof(double));

  /********************************************************************************/

  sprintf(Filename, "INPUT_STRUCTURE_%d", struct_index);
  Read_config_simple_3d(Num_particles, &fr_ref, Filename);

  sprintf(Filename, "minimized_struct_on_cusp");
  Read_config_simple_3d(Num_particles, &fr_ts, Filename);

  sprintf(Filename, "REFERENCE_FOR_ALIGNMENT");
  Read_config_simple_3d(Num_particles, &fr_ref_for_alignment, Filename);
  for(i = 0 ; i < Num_particles ; i++)
    {
      Reference_x[i+1] = fr_ref_for_alignment.x[i][0];
      Reference_y[i+1] = fr_ref_for_alignment.x[i][1];
      Reference_z[i+1] = fr_ref_for_alignment.x[i][2];
    }

  //Align_two_structures_lei(Reference_x, Reference_y, Reference_z, &fr_ref);
  //Align_two_structures_lei(Reference_x, Reference_y, Reference_z, &fr_ts);

  Calc_pair_distances_one_structure(Num_particles, fr_ref, pair_distances_structure);

  for(i = 0 ; i < Num_particles ; i++)
    {
      fr_start.x[i][0] = fr_ts.x[i][0];
      fr_start.x[i][1] = fr_ts.x[i][1];
      fr_start.x[i][2] = fr_ts.x[i][2];
    }


  is_step_size_adequate = 0;
  while(is_step_size_adequate == 0)
    {
       is_step_size_adequate = Check_step_size_offset(Num_particles, force_constant, cutoff, energy_offset, pair_distances_structure, Reference_x, Reference_y, Reference_z, rmsd_pathway, delta_rmsd_zero, &fr_start, step_size);
       if(is_step_size_adequate == 1)
         {
           printf("Step-size is %lf and it is adequate\n", step_size);
         }
       else
         {
           printf("Step-size is %lf and it is not adequate\n", step_size);
           step_size = step_size * step_size_reduction_factor;
         }
    }
  

  getforces_enm_offset(&fr_ref, pair_distances_structure, cutoff, force_constant, energy_offset);
  energy_ref = fr_ref.U;

  Num_struct_written = 0;
  Slide_down_single_step_size_offset(Num_particles, force_constant, cutoff, energy_offset, pair_distances_structure, Max_num_iterations, Reference_x, Reference_y, Reference_z, struct_index, rmsd_pathway, delta_rmsd_zero, &fr_ref, &fr_start, step_size, energy_ref, energy_from_ref_tol, &Num_struct_written, &fr_end, &iteration_number_end, step_size_reduction_factor);

  sprintf(Filename, "num_structures_written_%d", struct_index);
  out_num = fopen(Filename, "w");
  fprintf(out_num, "%d\n", Num_struct_written);
  fclose(out_num);

  getforces_enm_offset(&fr_end, pair_distances_structure, cutoff, force_constant, energy_offset);
  rmsd_from_ref_final = Calc_rmsd_one_structure(&fr_ref, &fr_end);
  sprintf(Filename, "final_struct_info_%d", struct_index);
  out = fopen(Filename, "w");
  fprintf(out, "Energy of the reference strucutre: %.8f\n", energy_ref);
  fprintf(out, "Number of iterations: %lu\n", iteration_number_end);
  fprintf(out, "Energy of final strucutre from minimization: %.8f\n", fr_end.U);
  fprintf(out, "RMSD of final strucutre from minimization from the reference strucutre: %.8f\n", rmsd_from_ref_final);
  fclose(out);
  sprintf(Filename, "final_struct_from_min_%d", struct_index);
  Write_config_simple(&fr_end, Filename); 

  return 0;
}












