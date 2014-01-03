/* Codes for calculation of rmsd of the string from initial string and/or string from the previous
 * iteration. 
 * Also provides codes for calculation of rmsd of structures*/

/* last modified on January 24, 2012*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "calc_rmsd.h"

double Calc_rmsd_one_structure(FR_DAT *fr, FR_DAT *fr_fixed_state)
{
  double msd, rmsd, DEL_X, DEL_Y, DEL_Z;
  int i;
  msd = 0.0;
  for(i = 0 ; i < fr->natoms; i++)
    {
      DEL_X = fr->x[i][0] - fr_fixed_state->x[i][0];
      DEL_Y = fr->x[i][1] - fr_fixed_state->x[i][1];
      DEL_Z = fr->x[i][2] - fr_fixed_state->x[i][2];
      msd = msd + pow(DEL_X, 2) + pow(DEL_Y, 2) + pow(DEL_Z, 2);
    }
  msd = msd / ((double)fr->natoms);
  rmsd = sqrt(msd);
  return rmsd;
}

/*Wrapper function.*/
void Calc_rmsd_all_structure(int Num_atoms, int Num_images, char *Filename_base, char *Filename_fixed_base, FILE *out_rmsd_structure, long int iteration_number)
{
  FR_DAT fr, fr_fixed_state;
  char Filename[1000];
  double rmsd;
  int i;

  fprintf(out_rmsd_structure, "%lu    ", iteration_number);
  for(i = 0 ; i < Num_images ; i++)
    {
      sprintf(Filename, "%s_%d", Filename_base, (i+1));
      Read_alloc_config_simple_3d(Num_atoms, &fr, Filename);
 
      sprintf(Filename, "%s_%d", Filename_fixed_base, (i+1));
      Read_alloc_config_simple_3d(Num_atoms, &fr_fixed_state, Filename);

      rmsd = Calc_rmsd_one_structure(&fr, &fr_fixed_state);
      if(i < (Num_images - 1))
        {
          fprintf(out_rmsd_structure, "%lf ", rmsd);
        }
      else
        {
          fprintf(out_rmsd_structure, "%lf\n", rmsd);
        }
    }

  /*Free memory*/
  free(fr.x);
  free(fr_fixed_state.x);
}


double Calc_rmsd_one_image_of_string(int Image_dimension, double *Image_position, double *Image_position_fixed)
{
  double msd, rmsd, DEL;
  int i;
  msd = 0.0;
  for(i = 0 ; i < Image_dimension ; i++)
    {
      DEL = Image_position[i] - Image_position_fixed[i];
      msd = msd + pow(DEL, 2);
    }
  msd = msd / ((double)Image_dimension);
  rmsd = sqrt(msd);
  return rmsd;
}

/* This function writes to two files. 
 * 1. Rmsd each image is stored.
 * 2. Average rmsd of the string is stored.
 * */
void Calc_rmsd_string(int Num_images, int Image_dimension, char *Filename_base, char *Filename_fixed_base, FILE *out_all_rmsd, FILE *out_average_rmsd, long int iteration_number)
{
  double *Image_position, *Image_position_fixed;
  char Filename[1000];
  double rmsd, average_rmsd;
  int i, d;
  FILE *in;

  Image_position = (double*)malloc(Image_dimension * sizeof(double));
  Image_position_fixed = (double*)malloc(Image_dimension * sizeof(double));  

  fprintf(out_all_rmsd, "%lu    ", iteration_number);
  fprintf(out_average_rmsd, "%lu    ", iteration_number);
  average_rmsd = 0.0;
  for(i = 0 ; i < Num_images ; i++)
    {
      sprintf(Filename, "%s_%d", Filename_base, (i+1));
      in = fopen(Filename, "r");
      if(in == NULL)
        {
          printf("File '%s' not found\n", Filename);
          exit(1);
        }
      for(d = 0 ; d < Image_dimension ; d++)
        {
          fscanf(in, "%lf ", &Image_position[d]);
        }
      fclose(in);

      sprintf(Filename, "%s_%d", Filename_fixed_base, (i+1));
      in = fopen(Filename, "r");
      if(in == NULL)
        {
          printf("File '%s' not found\n", Filename);
          exit(1);
        }
      for(d = 0 ; d < Image_dimension ; d++)
        {
          fscanf(in, "%lf ", &Image_position_fixed[d]);
        }
      fclose(in);

      rmsd = Calc_rmsd_one_image_of_string(Image_dimension, Image_position, Image_position_fixed);
      if(i < (Num_images - 1))
        {
          fprintf(out_all_rmsd, "%lf ", rmsd);
        }
      else
        {
          fprintf(out_all_rmsd, "%lf\n", rmsd);
        }

      average_rmsd = average_rmsd + rmsd;
    }
  average_rmsd = average_rmsd / ((double)Num_images);
  fprintf(out_average_rmsd, "%lf\n", average_rmsd); 
  fflush(out_average_rmsd); 

  free(Image_position);
  free(Image_position_fixed);
}


