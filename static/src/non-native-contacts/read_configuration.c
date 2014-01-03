/* Reads starting configuration from a file. The input file is a text file. This file only
   contains coordinate information. Velocities and forces are not stored.*/


/* last modified on January 19, 2012*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "read_configuration.h"


/* Format of the input file:
---------------------------------------------------------------
Number of atoms
Box length (one number, cubic box)
<blank line>
x_coordinate_of_particle_1 y_coordinate_of_particle_1 z_coordinate_of_particle_1
x_coordinate_of_particle_2 y_coordinate_of_particle_2 z_coordinate_of_particle_2
etc. total no. of lines should be equal to total number of atoms
*/
/* Reads the config file. It is assumed that required memory has already been allocated.
   Performs a consistency check to see whether the total numbers of atoms are consistent*/
void Read_config(FR_DAT *fr, char *input_filename, int Num_particles_from_top)
{
  int i, d;
  int dim = DIMENSION;
  FILE *in;
  double box_length_cubic;  

  in = fopen(input_filename, "r");
  if(in == NULL)
    {
      printf("File '%s' not found\n", input_filename);
      exit(1);
    }
  
  fscanf(in, "%d ", &fr->natoms);
  if(Num_particles_from_top != fr->natoms)
    {
      printf("Number of particles in input cofiguration file is different from that in topology file\n");
      exit(1);
    }

  fscanf(in, "%lf ", &box_length_cubic);
  for(i = 0 ; i < dim ; i++)
    {
      fr->box[i][i] = box_length_cubic;
      fr->box_dim2[i] = box_length_cubic / 2;
    }

  /* read particle coordinates*/
  for(i = 0 ; i < fr->natoms ; i++)
    {
      //fscanf(in, "%lf %lf %lf ", &fr->x[i][0], &fr->x[i][1], &fr->x[i][2]);
      for(d = 0 ; d < dim ; d++)
        {
          fscanf(in, "%lf ", &fr->x[i][d]);
        }
    }

  fclose(in);

}



/* Read the initial configuration and allocate necessary memory
 * Format of the initial configuration file is same as above.
 */
void Read_config_alloc(FR_DAT *fr, char *input_filename)
{
  int i, d;
  int dim = DIMENSION;
  FILE *in;
  double box_length_cubic;  

  in = fopen(input_filename, "r");
  if(in == NULL)
    {
      printf("File '%s' not found\n", input_filename);
      exit(1);
    }
  
  fscanf(in, "%d ", &fr->natoms);
  
  /*allocate memory*/
  fr->x = (rvec *)malloc(fr->natoms * sizeof(rvec));
  fr->f = (rvec *)malloc(fr->natoms * sizeof(rvec));
  fr->v = (rvec *)malloc(fr->natoms * sizeof(rvec));
  fr->m = (double *)malloc(fr->natoms * sizeof(double));

  fscanf(in, "%lf ", &box_length_cubic);
  for(i = 0 ; i < dim ; i++)
    {
      fr->box[i][i] = box_length_cubic;
      fr->box_dim2[i] = box_length_cubic / 2;
    }

  /* read particle coordinates*/
  for(i = 0 ; i < fr->natoms ; i++)
    {
      //fscanf(in, "%lf %lf %lf ", &fr->x[i][0], &fr->x[i][1], &fr->x[i][2]);
      for(d = 0 ; d < dim ; d++)
        {
          fscanf(in, "%lf ", &fr->x[i][d]);
        }
    }

  fclose(in);

}


/* Read the initial config file in simple format. No box length, the file consists of only 
 * coordinates of all the atoms/sites. The number of atoms read from the file is checked with
 * the number of atoms input from other source. This is suitable for non-periodic system 
 * where box length is irrelevant. It also allocates the required memory
 * Format:
 * Each line should have three numbers separated by white space giving x, y, and z cooridnates of one atom.
 * */
void Read_alloc_config_simple_3d(int Num_atoms, FR_DAT *fr, char *input_filename)
{
  //int d;
  //int dim = DIMENSION;
  FILE *in;
  char line[1000000];
  int Num_lines;
  //double val;
  double x, y, z;
  
  in = fopen(input_filename, "r");
  if(in == NULL)
    {
      printf("File '%s' not found\n", input_filename);
      exit(1);
    }

  fr->natoms = Num_atoms;
  fr->x = (rvec *)malloc(fr->natoms * sizeof(rvec));
  fr->f = (rvec *)malloc(fr->natoms * sizeof(rvec));
  fr->v = (rvec *)malloc(fr->natoms * sizeof(rvec));
  fr->m = (double *)malloc(fr->natoms * sizeof(double));

  /*Read coordinates*/
  Num_lines = 0;
  while(fgets(line, 1000000, in) != NULL)
    {
      if((Num_lines + 1) > Num_atoms)
        {
          printf("Number of atoms in file '%s' is more than the expected number of %d. Ignoring stuff after line %d\n", input_filename, Num_atoms, Num_lines);
          return;
        }
      /*for(d = 0 ; d < dim ; d++)
        {
          sscanf(line, "%lf", &val);
          fr->x[Num_lines][d] = val;
        }*/
      //sscanf(line, "%lf  %lf  %lf", &x, &y, &z);
      //one_coordinate_set->x[Num_lines][0] = x;
      //one_coordinate_set->x[Num_lines][1] = y;
      //one_coordinate_set->x[Num_lines][2] = z;
      sscanf(line, "%lf  %lf  %lf", &x, &y, &z);
      fr->x[Num_lines][0] = x;
      fr->x[Num_lines][1] = y;
      fr->x[Num_lines][2] = z;
      Num_lines = Num_lines + 1;
    }
  fclose(in); 
  /*Check if number of lines is less than number of atoms*/
  if(Num_lines < Num_atoms)
    {
      printf("Number of lines in file '%s' (%d) is less than the expected number of atoms (%d).\n", input_filename, Num_lines, Num_atoms);
      exit(1);
    }
  else
    {
      return;
    } 
}



/* Read the initial config file in simple format. No box length, the file consists of only 
 * coordinates of all the atoms/sites. The number of atoms read from the file is checked with
 * the number of atoms input from other source. This is suitable for non-periodic system 
 * where box length is irrelevant.
 * DOES NOT ALLOCATE REQUIRED MEMORY.
 * Format:
 * Each line should have three numbers separated by white space giving x, y, and z cooridnates of one atom.
 * THIS IS APPLICABLE FOR THREE DIMENSIONAL SYSTEMS.
 * */
void Read_config_simple_3d(int Num_atoms, FR_DAT *fr, char *input_filename)
{
  //int d;
  //int dim = DIMENSION;
  FILE *in;
  char line[1000000];
  int Num_lines;
  //double val;
  double x, y, z;
  
  in = fopen(input_filename, "r");
  if(in == NULL)
    {
      printf("File '%s' not found\n", input_filename);
      exit(1);
    }

  /*fr->natoms = Num_atoms;
  fr->x = (rvec *)malloc(fr->natoms * sizeof(rvec));
  fr->f = (rvec *)malloc(fr->natoms * sizeof(rvec));
  fr->v = (rvec *)malloc(fr->natoms * sizeof(rvec));
  fr->m = (double *)malloc(fr->natoms * sizeof(double));*/

  /*Read coordinates*/
  Num_lines = 0;
  while(fgets(line, 1000000, in) != NULL)
    {
      if((Num_lines + 1) > Num_atoms)
        {
          printf("Number of atoms in file '%s' is more than the expected number of %d. Ignoring stuff after line %d\n", input_filename, Num_atoms, Num_lines);
          return;
        }
      /*for(d = 0 ; d < dim ; d++)
        {
          sscanf(line, "%lf", &val);
          fr->x[Num_lines][d] = val;
        }*/
      sscanf(line, "%lf  %lf  %lf", &x, &y, &z);
      fr->x[Num_lines][0] = x;
      fr->x[Num_lines][1] = y;
      fr->x[Num_lines][2] = z;
      //one_coordinate_set->x[Num_lines][0] = x;
      //one_coordinate_set->x[Num_lines][1] = y;
      //one_coordinate_set->x[Num_lines][2] = z;
      Num_lines = Num_lines + 1;
    }
  fclose(in); 
  /*Check if number of lines is less than number of atoms*/
  if(Num_lines < Num_atoms)
    {
      printf("Number of lines in file '%s' (%d) is less than the expected number of atoms (%d).\n", input_filename, Num_lines, Num_atoms);
      exit(1);
    }
  else
    {
      return;
    } 
}



/* Read a VMD style xyz file.
 * Format:
 * First line: number of atom
 * Second line: a comment
 * Then for each line
 * atom name (C) x, y, z 
 * */
void Read_config_vmd_xyz(FR_DAT *fr, char *input_filename)
{
  FILE *in;
  int i;
  char line[1000000], str[100];

  in = fopen(input_filename, "r");
  if(in == NULL)
    {
      printf("File '%s' not found.\n", input_filename);
      exit(1);
    }
  fgets(line, 1000000, in);
  sscanf(line, "%d", &fr->natoms);
  fgets(line, 1000000, in);
  for(i = 0 ; i < fr->natoms ; i++)
    {
      fgets(line, 1000000, in);
      sscanf(line, "%s %lf %lf %lf", str, &fr->x[i][0], &fr->x[i][1], &fr->x[i][2]);
    }
  fclose(in);
}



/* Read a VMD style xyz file.
 * Format:
 * -------------------------
 * First line: number of atom
 * Second line: a comment
 * Then for each line
 * atom name (C) x, y, z 
 * -----------------------
 * Allocate memory.
 * */
void Read_alloc_config_vmd_xyz(FR_DAT *fr, char *input_filename)
{
  FILE *in;
  int i;
  char line[1000000], str[100];

  in = fopen(input_filename, "r");
  if(in == NULL)
    {
      printf("File '%s' not found.\n", input_filename);
      exit(1);
    }
  fgets(line, 1000000, in);
  sscanf(line, "%d", &fr->natoms);
  fr->x = (rvec*)malloc(fr->natoms * sizeof(rvec));
  fgets(line, 1000000, in);
  for(i = 0 ; i < fr->natoms ; i++)
    {
      fgets(line, 1000000, in);
      sscanf(line, "%s %lf %lf %lf", str, &fr->x[i][0], &fr->x[i][1], &fr->x[i][2]);
    }
  fclose(in);
}
