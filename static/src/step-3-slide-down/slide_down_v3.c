/* Routines for writing structures while sliding down a well.*/

/* Allows energy offsets*/

/* Last modified on May 20, 2013*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "slide_down_v3.h"



/* Check if the step-size if adequate. Return 1 if works, 0 if does not work*/
int Check_step_size_offset(int Num_particles, double force_constant, double cutoff, double energy_offset, double **pair_distances_structure, double *Reference_x, double *Reference_y, double *Reference_z, double rmsd_pathway, double delta_rmsd_zero, FR_DAT *fr_start, double step_size)
{
  FR_DAT fr;
  int i;
  double rmsd;
  int Num_structures_written;

  fr.natoms = Num_particles;
  fr.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr.f = (rvec*)malloc(Num_particles * sizeof(rvec));
  
  /* Copy the initial structure*/
  for(i = 0 ; i < Num_particles ; i++)
    {
      fr.x[i][0] = fr_start->x[i][0];
      fr.x[i][1] = fr_start->x[i][1];
      fr.x[i][2] = fr_start->x[i][2];
    }
  
  /* Calculate force*/
  getforces_enm_offset(&fr, pair_distances_structure, cutoff, force_constant, energy_offset);
  /* Move by one step*/
  for(i = 0 ; i < Num_particles ; i++)
    {
      fr.x[i][0] = fr.x[i][0] + (step_size * fr.f[i][0]);
      fr.x[i][1] = fr.x[i][1] + (step_size * fr.f[i][1]);
      fr.x[i][2] = fr.x[i][2] + (step_size * fr.f[i][2]);
    }
  Num_structures_written = 0;
  /* Check the rmsd of the present strucutre fro mthe initial strucutre. IF rmsd is greater than
     the target rmsd separation of the pathway then the step-size is clearly large*/
  Align_two_structures_lei(Reference_x, Reference_y, Reference_z, &fr);
  rmsd = Calc_rmsd_one_structure(&fr, fr_start);
  if(rmsd > rmsd_pathway) //rmsd is greated than the target, step-size inadequate
    {
      return 0;
    }
  else if(fabs(rmsd - rmsd_pathway) < delta_rmsd_zero) // write a strucutre, step-size adequate
    {
      Num_structures_written = Num_structures_written + 1;
      return 1;
    }
  else
    {
      while(rmsd < rmsd_pathway)
	{
	  getforces_enm_offset(&fr, pair_distances_structure, cutoff, force_constant, energy_offset);
	  for(i = 0 ; i < Num_particles ; i++)
	    {
	      fr.x[i][0] = fr.x[i][0] + (step_size * fr.f[i][0]);
	      fr.x[i][1] = fr.x[i][1] + (step_size * fr.f[i][1]);
	      fr.x[i][2] = fr.x[i][2] + (step_size * fr.f[i][2]);
	    }
	  Align_two_structures_lei(Reference_x, Reference_y, Reference_z, &fr);
	  rmsd = Calc_rmsd_one_structure(&fr, fr_start);
	  if(fabs(rmsd - rmsd_pathway) < delta_rmsd_zero)
	    {
	      Num_structures_written = Num_structures_written + 1;
	    }
	}
      /* RMSD of the current state is greater than the target rmsd, by this at least one strucutre should 
         have been written if the step-size was adequate. IF number of written structures is still zero then
         the step size is still not adequate.*/
      if(Num_structures_written == 0)
	{
	  return 0;
	}
      else
	{
	  return 1;
	}
    }

  free(fr.x);
  free(fr.f);
}


/* Slide down a well with a single step size.*/
void Slide_down_single_step_size_offset(int Num_particles, double force_constant, double cutoff, double energy_offset, double **pair_distances_structure, long int Max_num_iterations, double *Reference_x, double *Reference_y, double *Reference_z, int struct_index, double rmsd_pathway, double delta_rmsd_zero, FR_DAT *fr_ref, FR_DAT *fr_start, double step_size, double energy_ref, double energy_from_ref_tol, int *Num_struct_written, FR_DAT *fr_end, long int *iteration_number_end, double step_size_reduction_factor)
{
  FR_DAT fr, fr_last;
  int i;
  double rmsd_from_last_written, rmsd_from_ref;
  long int iteration_number;
  int is_complete;
  double energy_current;
  char Filename[1000];
  FILE *out_rmsd, *out_energy, *out;
  int is_step_size_adequate;

  fr.natoms = Num_particles;
  fr.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr.f = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_last.natoms = Num_particles;
  fr_last.x = (rvec*)malloc(Num_particles * sizeof(rvec));

  is_complete = 0;

  while(is_complete == 0)
    {
      printf("Trying step-size: %lf\n", step_size);
      is_step_size_adequate = 1;
      /* Copy the initial structure*/
      for(i = 0 ; i < Num_particles ; i++)
	{
	  fr.x[i][0] = fr_start->x[i][0];
	  fr.x[i][1] = fr_start->x[i][1];
	  fr.x[i][2] = fr_start->x[i][2];
	  
	  fr_last.x[i][0] = fr_start->x[i][0];
	  fr_last.x[i][1] = fr_start->x[i][1];
	  fr_last.x[i][2] = fr_start->x[i][2];
	}
      
      out = fopen("rmsd_all_steps", "w");
      sprintf(Filename, "rmsds_of_written_structures_%d", struct_index);
      out_rmsd = fopen(Filename, "w");
      sprintf(Filename, "energy_of_written_structures_%d", struct_index);
      out_energy = fopen(Filename, "w");
      iteration_number = 0;
      *Num_struct_written = 0;
      getforces_enm_offset(&fr, pair_distances_structure, cutoff, force_constant, energy_offset);
      energy_current = fr.U;
      do
	{ 
	  for(i = 0 ; i < Num_particles ; i++)
	    {
	      fr.x[i][0] = fr.x[i][0] + (step_size * fr.f[i][0]);
	      fr.x[i][1] = fr.x[i][1] + (step_size * fr.f[i][1]);
	      fr.x[i][2] = fr.x[i][2] + (step_size * fr.f[i][2]);
	    }
	  Align_two_structures_lei(Reference_x, Reference_y, Reference_z, &fr);
	  
	  rmsd_from_last_written = Calc_rmsd_one_structure(&fr, &fr_last);
	  rmsd_from_ref = Calc_rmsd_one_structure(&fr, fr_ref);
	  
	  fprintf(out, "%lf  %lf\n", rmsd_from_last_written, rmsd_from_ref);
	  
	  /* Check whether to write a strucutre file*/
	  if(fabs(rmsd_from_last_written - rmsd_pathway) < delta_rmsd_zero)
	    {
	      *Num_struct_written = *Num_struct_written + 1;
	      sprintf(Filename, "OUT_COORDS_SURFACE_%d_%d", struct_index, *Num_struct_written);
	      Write_config_simple(&fr, Filename);
	      fprintf(out_rmsd, "%d      %.8f    %.8f\n", *Num_struct_written, rmsd_from_last_written, rmsd_from_ref);
	      fflush(out_rmsd);
	      fprintf(out_energy, "%d   %lf\n", *Num_struct_written, fr.U);
	      fflush(out_energy);
	      for(i = 0 ; i < Num_particles ; i++)
		{
		  fr_last.x[i][0] = fr.x[i][0];
		  fr_last.x[i][1] = fr.x[i][1];
		  fr_last.x[i][2] = fr.x[i][2];
		}
	      
	    }
	  
	  if((rmsd_from_last_written > rmsd_pathway) && (fabs(rmsd_from_last_written - rmsd_pathway) > delta_rmsd_zero)) 
	    {
	      printf("RMSD from last written strucutre is greater than target rmsd for the pathway. Step-size is too large. %lf.\n", rmsd_from_last_written);
	      fclose(out_rmsd);
	      fclose(out_energy);
	      fclose(out);
	      //exit(1);
	      is_step_size_adequate = 0;
	      step_size = step_size * step_size_reduction_factor;
	    }
	  
	  
	  /* Figure out whether to stop the calculation*/
	  if((fabs(rmsd_from_ref - rmsd_pathway) < delta_rmsd_zero) || (fabs(energy_ref - energy_current) < energy_from_ref_tol))
	    {
	      is_complete = 1;
	    }
	  
	  /* Calculate force */
	  getforces_enm_offset(&fr, pair_distances_structure, cutoff, force_constant, energy_offset);
	  energy_current = fr.U;
	  
	  iteration_number = iteration_number + 1;
	}while((is_step_size_adequate == 1) && (is_complete == 0) && (iteration_number <= Max_num_iterations));
      
    }
  
  fclose(out_rmsd);
  fclose(out_energy);
  fclose(out);
  
  for(i = 0 ; i < Num_particles ; i++)
    {
      fr_end->x[i][0] = fr.x[i][0];
      fr_end->x[i][1] = fr.x[i][1];
      fr_end->x[i][2] = fr.x[i][2];
    }
  
  *iteration_number_end = iteration_number;
  
  free(fr.x);
  free(fr.f);
  free(fr_last.x);
}
