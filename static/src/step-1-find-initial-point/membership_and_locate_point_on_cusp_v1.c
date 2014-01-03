/* Locate point on the cusp hypersurface by iterative linear interpolation.*/

/* last update on May 20, 2013*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "membership_and_locate_point_on_cusp_v1.h"


/* Given a strucutre finds memebership on the two state harmonic potential with a cusp.
   Returns 1 if on surface 1.
   Returns 2 if on surface 2
   Returns 0 if on the cusp hypersurface, i.e. energies from two surfaces are same.
 */
int Find_membership_offset(FR_DAT *fr, double **pair_distances_structure_1, double cutoff_1, double force_constant_1, double energy_offset_1, double **pair_distances_structure_2, double cutoff_2, double force_constant_2, double energy_offset_2, double tol)
{
  double energy_1, energy_2;
  
  energy_1 = getforces_enm_return_pot_offset(fr, pair_distances_structure_1, cutoff_1, force_constant_1, energy_offset_1);
  energy_2 = getforces_enm_return_pot_offset(fr, pair_distances_structure_2, cutoff_2, force_constant_2, energy_offset_2);

  if(fabs(energy_1 - energy_2) < tol) /* Energies are same*/
    {
      return 0;
    }
  else if(energy_1 < energy_2)
    {
      return 1;
    }
  else
    {
      return 2;
    }
}




/* Generate strucutre from parameter. This is used for generating a new strucutre from the value 
   of the interpolating parameter. The end point structures used for interpolation are 
   'fr_1' and 'fr_2' respectively.
   Newly generated strucutre is stored in 'fr_new'.
 */
void Gen_struct_from_param(int Num_atoms, FR_DAT *fr_1, FR_DAT *fr_2, double t, FR_DAT *fr_new)
{
  double t_1, t_2;
  double r_1, r_2;
  int i, d;
  
  t_1 = 0.0;
  t_2 = 1.0;
  for(i = 0 ; i < Num_atoms ; i++)
    {
      for(d = 0 ; d < 3 ; d++)
	{
	  r_1 = fr_1->x[i][d];
	  r_2 = fr_2->x[i][d];
	  fr_new->x[i][d] = (r_1 * ((t-t_2) / (t_1-t_2))) + (r_2 * ((t-t_1) / (t_2-t_1)));
	}
    }
}


void Copy_config(int Num_atoms, FR_DAT *fr_old, FR_DAT *fr_new)
{
  int i, d;
  for(i = 0 ; i < Num_atoms ; i++)
    {
      for(d = 0 ; d < 3 ; d++)
        {
          fr_new->x[i][d] = fr_old->x[i][d];
        }
    }
}



/* Generate strucutre on the cusp. Search for strucutre on the cusp by recursively generating 
   two strucutres on two surfaces and checking the energy difference between them. When the 
   energy difference is smaller than the tolerance then two strucutres are said to be of equal 
   energy.
   End structures are stored in 'fr_1' and 'fr_2'.
   Num_images is the number of strucutres used for interpolation (including the end strucutres).
   'fr_try_1', 'fr_try_2' and 'fr_new' are used for calculation, these data strucutres 
   need to be allocated before passing on. This saves multiple memory allocation/deletion.*/
void Find_strucutre_on_cusp_offset(int Num_atoms, FR_DAT *fr_end_1, FR_DAT *fr_end_2, double tol, int Num_images, double **pair_distances_structure_1, double cutoff_1, double force_constant_1, double energy_offset_1, double **pair_distances_structure_2, double cutoff_2, double force_constant_2, double energy_offset_2, FR_DAT *fr_1, FR_DAT *fr_2, FR_DAT *fr_try_1, FR_DAT *fr_try_2, FR_DAT *fr_new, FR_DAT *fr_cusp)
{
  int alpha;
  int membership_previous, membership_current;
  int is_converged;
  double t;
  
  Copy_config(Num_atoms, fr_end_1, fr_1);
  Copy_config(Num_atoms, fr_end_2, fr_2);
  is_converged = 0;
  do
    {
      membership_previous = 1;
      for(alpha = 1 ; alpha < (Num_images-1) ; alpha++)
        {
          t = ((double)alpha) / ((double)(Num_images - 1));
	  Gen_struct_from_param(Num_atoms, fr_1, fr_2, t, fr_new);
	  membership_current = Find_membership_offset(fr_new, pair_distances_structure_1, cutoff_1, force_constant_1, energy_offset_1, pair_distances_structure_2, cutoff_2, force_constant_2, energy_offset_2, tol);
          if(membership_current == 0)   /* On the cusp*/
            {
	      Copy_config(Num_atoms, fr_new, fr_cusp);
              is_converged = 1;
            }
          else if(membership_current == membership_previous)
            {
              Copy_config(Num_atoms, fr_new, fr_try_1);
              continue;
            }
          else
            {
	      Copy_config(Num_atoms, fr_new, fr_try_2);
              break;
            }
        }
      Copy_config(Num_atoms, fr_try_1, fr_1);
      Copy_config(Num_atoms, fr_try_2, fr_2);
    }while(is_converged != 1);
}



