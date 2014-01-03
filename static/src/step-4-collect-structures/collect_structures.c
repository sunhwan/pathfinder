/* Collects structures for constructing the pathway.*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "collect_structures.h"

void Collect_structures(int Num_particles, int Num_structures_1, int Num_structures_2, int *Num_struct_in_pathway)
{
  FR_DAT fr;
  int c, k;
  char Filename[1000];
  int index_of_ts_1, index_of_ts_2;
  FILE *out;

  fr.natoms = Num_particles;
  fr.x = (rvec*)malloc(Num_particles * sizeof(rvec));

  c = 0;
  /* Read all the structures*/
  /* We update the image index as we read.*/
  /* Read reference structure 1, This is image 0*/
  sprintf(Filename, "INPUT_STRUCTURE_1");
  Read_config_simple_3d(Num_particles, &fr, Filename);
  c = c + 1;
  sprintf(Filename, "COORDS_IMAGE_%d", c);
  Write_config_simple(&fr, Filename);

  /* Read all the structure from surface 1 in backwards the last created structure is read first*/
  for(k = 0 ; k < Num_structures_1 ; k++)
    {
      sprintf(Filename, "OUT_COORDS_SURFACE_1_%d", (Num_structures_1 - k));
      Read_config_simple_3d(Num_particles, &fr, Filename);
      c = c + 1;
      sprintf(Filename, "COORDS_IMAGE_%d", c);
      Write_config_simple(&fr, Filename);
    }
  /* Read the transition state on surface 1*/
  sprintf(Filename, "TS_1");
  Read_config_simple_3d(Num_particles, &fr, Filename);
  c = c + 1;
  index_of_ts_1 = c;
  sprintf(Filename, "COORDS_IMAGE_%d", c);
  Write_config_simple(&fr, Filename);
  /* Read the transition state on surface 2*/
  sprintf(Filename, "TS_2");
  Read_config_simple_3d(Num_particles, &fr, Filename);
  c = c + 1;
  index_of_ts_2 = c;
  sprintf(Filename, "COORDS_IMAGE_%d", c);
  Write_config_simple(&fr, Filename);
  /* Read all structures from surface 2*/
  for(k = 0 ; k < Num_structures_2 ; k++)
    {
      sprintf(Filename, "OUT_COORDS_SURFACE_2_%d", (k + 1));
      Read_config_simple_3d(Num_particles, &fr, Filename);
      c = c + 1;
      sprintf(Filename, "COORDS_IMAGE_%d", c);
      Write_config_simple(&fr, Filename);
    }
  /* Read the second reference structure. This is the last image i.e. (Num_images-1)*/
  sprintf(Filename, "INPUT_STRUCTURE_2");
  Read_config_simple_3d(Num_particles, &fr, Filename);
  c = c + 1;
  sprintf(Filename, "COORDS_IMAGE_%d", c);
  Write_config_simple(&fr, Filename);

  *Num_struct_in_pathway = c;

  out = fopen("number_of_structures_in_pathway", "w");
  fprintf(out, "Total number of structures in the pathway: %d\n", c);
  fprintf(out, "Index of TS 1: %d\n", index_of_ts_1);
  fprintf(out, "Index of TS 2: %d\n", index_of_ts_2);
  fclose(out);

  free(fr.x);
}



/* Collect strucutres when one transition state strucutre is used.*/
void Collect_structures_one_ts(int Num_particles, FR_DAT *fr_ts, int Num_structures_1, int Num_structures_2, int *Num_struct_in_pathway, int *index_of_ts_struct)
{
  FR_DAT fr;
  int c, k;
  char Filename[1000];
  FILE *out;

  fr.natoms = Num_particles;
  fr.x = (rvec*)malloc(Num_particles * sizeof(rvec));

  c = 0;
  /* Read all the structures*/
  /* We update the image index as we read.*/
  /* Read reference structure 1, This is image 0*/
  sprintf(Filename, "INPUT_STRUCTURE_1");
  Read_config_simple_3d(Num_particles, &fr, Filename);
  c = c + 1;
  sprintf(Filename, "COORDS_IMAGE_%d", c);
  Write_config_simple(&fr, Filename);

  /* Read all the structure from surface 1 in backwards the last created structure is read first*/
  for(k = 0 ; k < Num_structures_1 ; k++)
    {
      sprintf(Filename, "OUT_COORDS_SURFACE_1_%d", (Num_structures_1 - k));
      Read_config_simple_3d(Num_particles, &fr, Filename);
      c = c + 1;
      sprintf(Filename, "COORDS_IMAGE_%d", c);
      Write_config_simple(&fr, Filename);
    }
  /* Write the transition state structure*/
  c = c + 1;
  *index_of_ts_struct = c;
  sprintf(Filename, "COORDS_IMAGE_%d", c);
  Write_config_simple(fr_ts, Filename);
  /* Read all structures from surface 2*/
  for(k = 0 ; k < Num_structures_2 ; k++)
    {
      sprintf(Filename, "OUT_COORDS_SURFACE_2_%d", (k + 1));
      Read_config_simple_3d(Num_particles, &fr, Filename);
      c = c + 1;
      sprintf(Filename, "COORDS_IMAGE_%d", c);
      Write_config_simple(&fr, Filename);
    }
  /* Read the second reference structure. This is the last image i.e. (Num_images-1)*/
  sprintf(Filename, "INPUT_STRUCTURE_2");
  Read_config_simple_3d(Num_particles, &fr, Filename);
  c = c + 1;
  sprintf(Filename, "COORDS_IMAGE_%d", c);
  Write_config_simple(&fr, Filename);

  *Num_struct_in_pathway = c;

  out = fopen("number_of_structures_in_pathway", "w");
  fprintf(out, "Total number of structures in the pathway: %d\n", c);
  fprintf(out, "Index of TS: %d\n", *index_of_ts_struct);
  fclose(out);

  free(fr.x);
}


