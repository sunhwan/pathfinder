/* Starting from a strucutre on the cusp surface minimize it on the cusp surface.*/

/* Allows energy offsets for both the end states*/

/* Last modified on May 20, 2013*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "include_all.h"




/* Input file: 'INPUT_MINIMIZE_ON_CUSP'
 * Format:
 * -----------------------------------------
 *  NUM_PARTICLES (Number of particles in the system %d)
 *  CUT_OFFS (Cutoffs for both the end states. %lf %lf)
 *  FORCE_CONSTANTS (Force constants for both the end states. %lf %lf)
 *  ENERGY_OFFSETS (Energy offsets for both the end states. %lf %lf)
 *  CUSP_TOLERANCE (Energy tolerance for determining strucutre on cusp. %lf)
 *  NUM_IMAGES_INTERPOLATE (Number of images used for interpolation. %d)
 *  STEP_SIZES (Step sizes for minimization for two surfaces. %lf %lf)
 *  NUM_ITER_CHECK_STEP_SIZES (Number of iterations of initial run for checking step-sizes. %d)
 *  STEP_SIZE_REDUCTION_FACTORS (Factors by which step-sizes are decreased while searching the optimum step-sizes %lf %lf)
 *  ENERGY_TOL_CONVERGENCE (Covergence criterion based on energy. %lf) 
 *  NUM_ITERATIONS (Number of iterations. %lu)
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
  FR_DAT fr_end_1, fr_end_2, fr_1, fr_2, fr_try_1, fr_try_2, fr_new, fr_cusp, fr_cusp_ini;
  FILE *in, *out, *out_rmsd;
  int i;
  double potential_1, potential_2;
  double step_size_1, step_size_2;
  long int Num_iterations, iteration_number;
  FR_DAT fr_ref_for_alignment;
  double *Reference_x, *Reference_y, *Reference_z;
  double rmsd_from_ini, rmsd_cusp_end_1, rmsd_cusp_end_2;
  double energy_tol_convergence;
  double energy_previous_step;
  int is_converged;
  int Num_points;
  double step_size_reduction_factor_1, step_size_reduction_factor_2;

  in = fopen("INPUT_MINIMIZE_ON_CUSP", "r");
  if(in == NULL)
    {
      printf("File 'INPUT_MINIMIZE_ON_CUSP' not found.\n");
      exit(1);
    }
  fscanf(in, "%s %d ", key, &Num_particles);
  fscanf(in, "%s %lf %lf ", key, &cutoff_1, &cutoff_2);
  fscanf(in, "%s %lf %lf ", key, &force_constant_1, &force_constant_2);
  fscanf(in, "%s %lf %lf ", key, &energy_offset_1, &energy_offset_2);
  fscanf(in, "%s %lf ", key, &tol);
  fscanf(in, "%s %d ", key, &Num_images);
  fscanf(in, "%s %lf %lf ", key, &step_size_1, &step_size_2);
  fscanf(in, "%s %d ", key, &Num_points);
  fscanf(in, "%s %lf %lf ", key, &step_size_reduction_factor_1, &step_size_reduction_factor_2);
  fscanf(in, "%s %lf ", key, &energy_tol_convergence);
  fscanf(in, "%s %lu ", key, &Num_iterations); 
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
  fr_end_1.f = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_end_2.natoms = Num_particles ;
  fr_end_2.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_end_2.f = (rvec*)malloc(Num_particles * sizeof(rvec));
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
  fr_ref_for_alignment.natoms = Num_particles;
  fr_ref_for_alignment.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_cusp_ini.natoms = Num_particles;
  fr_cusp_ini.x = (rvec*)malloc(Num_particles * sizeof(rvec));

  Reference_x = (double*)malloc((Num_particles + 1) * sizeof(double));
  Reference_y = (double*)malloc((Num_particles + 1) * sizeof(double));
  Reference_z = (double*)malloc((Num_particles + 1) * sizeof(double));

  sprintf(Filename, "INPUT_STRUCTURE_1");
  Read_config_simple_3d(Num_particles, &fr_end_1, Filename);
  sprintf(Filename, "INPUT_STRUCTURE_2");
  Read_config_simple_3d(Num_particles, &fr_end_2, Filename);
 
  Calc_pair_distances_both_structures(Num_particles, fr_end_1, fr_end_2, pair_distances_structure_1, pair_distances_structure_2);

  sprintf(Filename, "REFERENCE_FOR_ALIGNMENT");
  Read_config_simple_3d(Num_particles, &fr_ref_for_alignment, Filename);
  for(i = 0 ; i < Num_particles ; i++)
    {
      Reference_x[i+1] = fr_ref_for_alignment.x[i][0];
      Reference_y[i+1] = fr_ref_for_alignment.x[i][1];
      Reference_z[i+1] = fr_ref_for_alignment.x[i][2];
    }


  /*******************************************************************************/

  /* Check and modify step-sizes*/
  Choose_step_size_offset(Num_particles, force_constant_1, force_constant_2, cutoff_1, cutoff_2, energy_offset_1, energy_offset_2, pair_distances_structure_1, pair_distances_structure_2, tol, Num_images, &fr_end_1, &fr_end_2, &fr_1, &fr_2, &fr_try_1, &fr_try_2, &fr_new, &fr_cusp, Reference_x, Reference_y, Reference_z, Num_points, step_size_reduction_factor_1, step_size_reduction_factor_2, &step_size_1, &step_size_2);

  //exit(1);

  /* Minimize on cusp*/
  sprintf(Filename, "initial_struct_on_cusp");
  Read_config_simple_3d(Num_particles, &fr_cusp, Filename);
  Read_config_simple_3d(Num_particles, &fr_cusp_ini, Filename);

  out = fopen("energies_along_minimization", "w");
  potential_1 = getforces_enm_return_pot_offset(&fr_cusp, pair_distances_structure_1, cutoff_1, force_constant_1, energy_offset_1);
  potential_2 = getforces_enm_return_pot_offset(&fr_cusp, pair_distances_structure_2, cutoff_2, force_constant_2, energy_offset_2);
  fprintf(out, "%d  %.15e  %.15e\n", 0, potential_1, potential_2);
  fflush(out);
  energy_previous_step = potential_1;
  out_rmsd = fopen("rmsd_from_initial", "w");
  fprintf(out_rmsd, "%d    %.8f\n", 0, 0.0);
  iteration_number = 1;
  is_converged = 0;
  do
    {
      /* Minimize for one step on both the surfaces*/
      /** Energy Surface 1**/
      Copy_config(Num_particles, &fr_cusp, &fr_end_1);
      getforces_enm_offset(&fr_end_1, pair_distances_structure_1, cutoff_1, force_constant_1, energy_offset_1);
      /** Energy Surface 2**/
      Copy_config(Num_particles, &fr_cusp, &fr_end_2);
      getforces_enm_offset(&fr_end_2, pair_distances_structure_2, cutoff_2, force_constant_2, energy_offset_2);
      /** Minimize**/
      for(i = 0 ; i < Num_particles ; i++)
        {
          fr_end_1.x[i][0] = fr_end_1.x[i][0] + (step_size_1 * fr_end_1.f[i][0]);
          fr_end_1.x[i][1] = fr_end_1.x[i][1] + (step_size_1 * fr_end_1.f[i][1]);
          fr_end_1.x[i][2] = fr_end_1.x[i][2] + (step_size_1 * fr_end_1.f[i][2]);

          fr_end_2.x[i][0] = fr_end_2.x[i][0] + (step_size_2 * fr_end_2.f[i][0]);
          fr_end_2.x[i][1] = fr_end_2.x[i][1] + (step_size_2 * fr_end_2.f[i][1]);
          fr_end_2.x[i][2] = fr_end_2.x[i][2] + (step_size_2 * fr_end_2.f[i][2]);
        }
 
      /* Align*/
      Align_two_structures_lei(Reference_x, Reference_y, Reference_z, &fr_end_1);
      Align_two_structures_lei(Reference_x, Reference_y, Reference_z, &fr_end_2);

      /* Create a point on the cusp by interpolating between two new points*/
      Find_strucutre_on_cusp_offset(Num_particles, &fr_end_1, &fr_end_2, tol, Num_images, pair_distances_structure_1, cutoff_1, force_constant_1, energy_offset_1, pair_distances_structure_2, cutoff_2, force_constant_2, energy_offset_2, &fr_1, &fr_2, &fr_try_1, &fr_try_2, &fr_new, &fr_cusp);

      /* Energy at the new point on the cusp*/
      potential_1 = getforces_enm_return_pot_offset(&fr_cusp, pair_distances_structure_1, cutoff_1, force_constant_1, energy_offset_1);
      potential_2 = getforces_enm_return_pot_offset(&fr_cusp, pair_distances_structure_2, cutoff_2, force_constant_2, energy_offset_2);
      fprintf(out, "%lu  %.15e  %.15e\n", iteration_number, potential_1, potential_2);    
      fflush(out);
      /* Calculate rmsd of the current point on the cusp from it's initial position*/
      rmsd_from_ini = Calc_rmsd_one_structure(&fr_cusp, &fr_cusp_ini);
      fprintf(out_rmsd, "%lu    %.8f\n", iteration_number, rmsd_from_ini);
      fflush(out_rmsd);
      //printf("%lu\n", iteration_number);
 
      /* Determine convergence*/
      if(fabs(potential_1 - energy_previous_step) < energy_tol_convergence)
        {
          is_converged = 1;
        }

      /* Reset energy for the previous set for checking convergence in the next step*/
      energy_previous_step = potential_1;
      
      /* Increase pointer*/
      iteration_number = iteration_number + 1;
    }while((is_converged == 0) && (iteration_number <= Num_iterations)); 

  fclose(out);
  fclose(out_rmsd);
  
  /* Write the point on cusp*/
  sprintf(Filename, "minimized_struct_on_cusp");
  Write_config_simple(&fr_cusp, Filename); 

  /* RMSDs of the TS fro mthe end structures*/
  sprintf(Filename, "INPUT_STRUCTURE_1");
  Read_config_simple_3d(Num_particles, &fr_end_1, Filename);
  sprintf(Filename, "INPUT_STRUCTURE_2");
  Read_config_simple_3d(Num_particles, &fr_end_2, Filename);
  rmsd_cusp_end_1 = Calc_rmsd_one_structure(&fr_cusp, &fr_end_1);
  rmsd_cusp_end_2 = Calc_rmsd_one_structure(&fr_cusp, &fr_end_2);
  out = fopen("minimization_on_cusp.log", "w");
  fprintf(out, "Final step-sizes: %lf   %lf\n", step_size_1, step_size_2);
  fprintf(out, "RMSD of the energy minimized structure on the cusp from end point 1: %lf\n", rmsd_cusp_end_1);
  fprintf(out, "RMSD of the energy minimized structure on the cusp from end point 2: %lf\n", rmsd_cusp_end_2);
  fclose(out); 

  return 0;
}



