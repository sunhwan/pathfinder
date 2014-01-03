/* Prepares an initial string from two end structures by linear interpolation.
 * This is for the case when all the coordinates of the system are used 
 * for representing the string.*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "prepare_initial_string_all_atom.h"


/* THIS IS THE CODE FOR A TWO DIMENSIONAL SYSTEM.
 * Initialize the string by distribution images at equal distances between the first and the last points 
 * along a straight line. Initial string is a straight line. */
void Init_string(FR_DAT *string)
{
  int Num_images = string->natoms;

  /*
   * Linear interpolation.
   * If for a parameter t
   * x(t_1) = x_1  y(t_1) = y_1 
   * x(t_2) = x_2  y(t_2) = y_2
   * Then the interpolation formulas that satisfies these conditions
   * x(t) = x_1 * [(t-t_2) / (t_1-t_2)] + x_2 * [(t-t_1) / (t_2-t_1)]
   * y(t) = y_1 * [(t-t_2) / (t_1-t_2)] + y_2 * [(t-t_1) / (t_2-t_1)]
   **/

   int i;
   double x_1 = string->x[0][0];
   double y_1 = string->x[0][1];
   double x_2 = string->x[Num_images -1 ][0];
   double y_2 = string->x[Num_images -1 ][1];
   double t_1 = 0.0;
   double t_2 = 1.0;
   double t;
   /*here t = 0 for (x_1, y_1) and t = 1 for (x_2,y_2). we want points equally spaced in t*/
   for(i = 1 ; i < (Num_images - 1) ; i++)  /*frist and last points are fixed*/
     {
       t = ((double)i) / ((double)(Num_images - 1));
       string->x[i][0] = x_1 * ((t-t_2) / (t_1-t_2)) + x_2 * ((t-t_1) / (t_2-t_1));
       string->x[i][1] = y_1 * ((t-t_2) / (t_1-t_2)) + y_2 * ((t-t_1) / (t_2-t_1));
     }
}


/* Given two end structures prepares all structures for all the images and write them to separate files.*/
void Prepare_initial_string_all_atom(int Num_atoms, char Filename_input_structure_1[1000], char Filename_input_structure_2[1000], int Num_images)
{
  int i, d, im;
  double t, t_1, t_2;
  double r_1, r_2;
  FR_DAT fr_1, fr_2, fr;
  char Filename_output[1000];

  Read_alloc_config_simple_3d(Num_atoms, &fr_1, Filename_input_structure_1);
  Read_alloc_config_simple_3d(Num_atoms, &fr_2, Filename_input_structure_2);
  fr.natoms = Num_atoms;
  fr.x = (rvec*)malloc(Num_atoms * sizeof(rvec));

  sprintf(Filename_output, "COORDS_IMAGE_1");  
  Write_config_simple(&fr_1, Filename_output);

  sprintf(Filename_output, "COORDS_IMAGE_%d", Num_images);
  Write_config_simple(&fr_2, Filename_output);

  t_1 = 0.0;
  t_2 = 1.0;

  FILE *out = fopen("VALUES_OF_PARAMETERS", "w");
  fprintf(out, "%d  %.15e\n", 1, 0.0);

  for(im = 1 ; im < (Num_images - 1) ; im++)
    {
      t = ((double)im) / ((double)(Num_images - 1)); 
      /*Loop over all the particles*/
      for(i = 0 ; i < Num_atoms ; i++)
        {
          for(d = 0 ; d < 3 ; d++)
            {
              r_1 = fr_1.x[i][d];
              r_2 = fr_2.x[i][d];
              fr.x[i][d] = (r_1 * ((t-t_2) / (t_1-t_2))) + (r_2 * ((t-t_1) / (t_2-t_1))); 
            }
        }
      sprintf(Filename_output, "COORDS_IMAGE_%d", (im + 1)); 
      Write_config_simple(&fr, Filename_output);
      fprintf(out, "%d  %.15e\n", (im + 1), t);
    }

  fprintf(out, "%d  %.15e\n", Num_images, 1.0);
  fclose(out);
  
}


