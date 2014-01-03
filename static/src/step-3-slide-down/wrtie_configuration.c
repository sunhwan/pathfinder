/* Writes a configuration file. Used for writing the final configuration at the end of the
   run. */

/* Last modified on November 01, 2011*/


/* Format of the output file:
---------------------------------------------------------------
Number of atoms
Box length (one number, cubic box)
<blank line>
x_coordinate_of_particle_1 y_coordinate_of_particle_1 z_coordinate_of_particle_1
x_coordinate_of_particle_2 y_coordinate_of_particle_2 z_coordinate_of_particle_2
etc. total no. of lines should be equal to total number of atoms
*/
/* Writes a config file. Only coordinates no velocities or forces*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "wrtie_configuration.h"

void Write_config(FR_DAT *fr, char *output_filename)
{
  int i, d;
  int dim = DIMENSION;
  FILE *out;
  double box_length_cubic;

  box_length_cubic = fr->box[0][0];

  out = fopen(output_filename, "w");

  fprintf(out, "%d\n%.15e\n\n", fr->natoms, box_length_cubic);

  /* write particle coordinates*/
  for(i = 0 ; i < fr->natoms ; i++)
    {
      //fprintf(out, "%.15e %.15e %.15e\n", fr->x[i][0], fr->x[i][1], fr->x[i][2]);
      for(d = 0 ; d < dim ; d ++)
        {
          fprintf(out, "%.15e ", fr->x[i][d]);
        }
      fprintf(out, "\n");
    }

  fclose(out);

}



/* Writes configurtions (frames) in the trajectory file. This file is used for post 
   processing. 

   The format of the file:
   ----------------------------------------------------------------------------------------
   No. of atoms
   No. of frames
   Box length
   <Blank line>
   Coordinates of forces for all the particles. One line for each particle with the follwoing format:
   coordinate_x   coordinate_y   coordinate_z    force_x  force_y  force_z
   ...  ...  ...
   <Blank line>
   same for the next frame
   ----------------------------------------------------------------------------------------

   This functions writes the infromation for one frame. Data at the begining of the file 
   is writen when the file is first created.
*/

void Write_config_to_trajectory_file(FR_DAT *fr, FILE *output_file_pointer)
{
  int i, d;
  int dim = DIMENSION;

  for(i = 0 ; i < fr->natoms ; i++)
    {
      //fprintf(output_file_pointer, "%.15e   %.15e   %.15e    %.15e   %.15e   %.15e\n", fr->x[i][0], fr->x[i][1], fr->x[i][2], fr->f[i][0], fr->f[i][1], fr->f[i][2]);
      for(d = 0 ; d < dim ; d ++)
        {
          fprintf(output_file_pointer, "%.15e ", fr->x[i][d]);
        }
      fprintf(output_file_pointer," ");
      for(d = 0 ; d < dim ; d ++)
        {
          fprintf(output_file_pointer, "%.15e ", fr->f[i][d]);
        }
      fprintf(output_file_pointer, "\n");
    }
  fprintf(output_file_pointer, "\n");
}



/* New format for trajectory file. forces are not stored, only coordiates and nothing else. 
 * This will help us to process the trajectory of a one particle model system*/
void Write_config_to_trajectory_file_only_coordinates(FR_DAT *fr, FILE *output_file_pointer)
{
  int i, d;
  int dim = DIMENSION;

  for(i = 0 ; i < fr->natoms ; i++)
    {
      //fprintf(output_file_pointer, "%.15e   %.15e   %.15e    %.15e   %.15e   %.15e\n", fr->x[i][0], fr->x[i][1], fr->x[i][2], fr->f[i][0], fr->f[i][1], fr->f[i][2]);
      for(d = 0 ; d < dim ; d ++)
        {
          fprintf(output_file_pointer, "%.15e ", fr->x[i][d]);
        }
      //fprintf(output_file_pointer," ");
      /*for(d = 0 ; d < dim ; d ++)
        {
          fprintf(output_file_pointer, "%.15e ", fr->f[i][d]);
	}
      */
      fprintf(output_file_pointer, "\n");
    }
  //fprintf(output_file_pointer, "\n");
}

/* Write config file with only coordinates.
 * Format:
 * Each line has the coordinates of an atom separated by white space, starting from x, y, ..etc.*/
void Write_config_simple(FR_DAT *fr, char *output_filename)
{
  int i, d;
  int dim = DIMENSION;
  FILE *out;
  
  out = fopen(output_filename, "w");
  
  /* write particle coordinates*/
  for(i = 0 ; i < fr->natoms ; i++)
    {
      for(d = 0 ; d < dim ; d++)
	{
	  fprintf(out, "%.15e ", fr->x[i][d]);
	}
      fprintf(out, "\n");
    }
  
  fclose(out);
  
}


/* Wrties a VMD style xyz file. This is needed for visualization.
 * Format:
 * First line: number of atom
 * Second line: a comment
 * Then for each line
 * atom name (C) x, y, z 
 * */
void Write_config_vmd_xyz(FR_DAT *fr, char *output_filename)
{
  int i;
  FILE *out;

  out = fopen(output_filename, "w");
  fprintf(out, "%d\n", fr->natoms);
  fprintf(out, " file generated by my program\n");
  /* write particle coordinates*/
  for(i = 0 ; i < fr->natoms ; i++)
    {
      fprintf(out, " C        %lf       %lf       %lf\n", fr->x[i][0], fr->x[i][1], fr->x[i][2]);
    }

  fclose(out);
}