/* Routines for checking the input step-sizes*/

/* Last modified on May 20, 2013*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "check_step_sizes.h"


/* Ckecks if the values of the energies are decreasing.
   Return 1 if values are decreasing.
   Return 0 if vlaues are not increasing*/
int Is_energy_decreasing(int Num_points, double *energy_values)
{
  int i;
  double previous_value;
  
  previous_value = energy_values[0];
  for(i = 1 ; i < Num_points ; i++)
    {
      if(energy_values[i] > previous_value)
       {
         return 0;
       }
      else
       {
         previous_value = energy_values[i];
       }
    }
  return 1;
}





/* Choose step-sizes by monitoring the energ values of a short run.
   If they are monotonically decreasing the step-sizes are adequate,
   otherwise not.*/
void Choose_step_size_offset(int Num_particles, double force_constant_1, double force_constant_2, double cutoff_1, double cutoff_2, double energy_offset_1, double energy_offset_2, double **pair_distances_structure_1, double **pair_distances_structure_2, double tol, int Num_images, FR_DAT *fr_end_1, FR_DAT *fr_end_2, FR_DAT *fr_1, FR_DAT *fr_2, FR_DAT *fr_try_1, FR_DAT *fr_try_2, FR_DAT *fr_new, FR_DAT *fr_cusp, double *Reference_x, double *Reference_y, double *Reference_z, int Num_points, double step_size_reduction_factor_1, double step_size_reduction_factor_2, double *step_size_1, double *step_size_2)
{
  int i, n;
  char Filename[1000];
  /*1, if they are 0 otherwise*/
  int Are_step_sizes_adequate;
  /*1, yes; 0, no*/
  int Is_energy_increasing;
  double potential_1, potential_previous;

  Are_step_sizes_adequate = 0;
  while(Are_step_sizes_adequate == 0)
    {
       Is_energy_increasing = 0;
       /* Perform minimization for few steps and collect energies*/
       sprintf(Filename, "initial_struct_on_cusp");
       Read_config_simple_3d(Num_particles, fr_cusp, Filename);
       potential_1 = getforces_enm_return_pot_offset(fr_cusp, pair_distances_structure_1, cutoff_1, force_constant_1, energy_offset_1);
       potential_previous = potential_1;
       printf("%lf\n", potential_1);
       for(n = 0 ; n < Num_points ; n++)
         {
           Copy_config(Num_particles, fr_cusp, fr_end_1);
           getforces_enm_offset(fr_end_1, pair_distances_structure_1, cutoff_1, force_constant_1, energy_offset_1); 
           Copy_config(Num_particles, fr_cusp, fr_end_2);
           getforces_enm_offset(fr_end_2, pair_distances_structure_2, cutoff_2, force_constant_2, energy_offset_2);
           for(i = 0 ; i < Num_particles ; i++)
             {
               fr_end_1->x[i][0] = fr_end_1->x[i][0] + ((*step_size_1) * fr_end_1->f[i][0]);
               fr_end_1->x[i][1] = fr_end_1->x[i][1] + ((*step_size_1) * fr_end_1->f[i][1]);
               fr_end_1->x[i][2] = fr_end_1->x[i][2] + ((*step_size_1) * fr_end_1->f[i][2]);

               fr_end_2->x[i][0] = fr_end_2->x[i][0] + ((*step_size_2) * fr_end_2->f[i][0]);
               fr_end_2->x[i][1] = fr_end_2->x[i][1] + ((*step_size_2) * fr_end_2->f[i][1]);
               fr_end_2->x[i][2] = fr_end_2->x[i][2] + ((*step_size_2) * fr_end_2->f[i][2]);
             }
             Align_two_structures_lei(Reference_x, Reference_y, Reference_z, fr_end_1);
             Align_two_structures_lei(Reference_x, Reference_y, Reference_z, fr_end_2);
             Find_strucutre_on_cusp_offset(Num_particles, fr_end_1, fr_end_2, tol, Num_images, pair_distances_structure_1, cutoff_1, force_constant_1, energy_offset_1, pair_distances_structure_2, cutoff_2, force_constant_2, energy_offset_2, fr_1, fr_2, fr_try_1, fr_try_2, fr_new, fr_cusp);
             potential_1 = getforces_enm_return_pot_offset(fr_cusp, pair_distances_structure_1, cutoff_1, force_constant_1, energy_offset_1);
             printf("%lf\n", potential_1);

             if(potential_1 > potential_previous)
               {
                 Is_energy_increasing = 1;
                 break;
               }
             else
               {
                 potential_previous = potential_1;
               }
         }


       printf("%d\n", n);
       if(Is_energy_increasing == 0)
         {
           printf("Step-sizes are adequate\n");
           Are_step_sizes_adequate = 1;
         }
       else
         {
           printf("Step-sizes are not adequate\n");
           *step_size_1 = (*step_size_1) * step_size_reduction_factor_1;
           *step_size_2 = (*step_size_2) * step_size_reduction_factor_2;
           printf("Trying new step-sizes: %lf  %lf\n", (*step_size_1), (*step_size_2));
         }
  
    }

}
