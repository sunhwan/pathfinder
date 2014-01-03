/* Writes the coordinates of images to a PDB trajectory file
 * The cooridnates should be in files with names 'COORDS_IMAGE_1', 'COORDS_IMAGE_2' etc.*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "write_images_to_PDB_traj_v2.h"


void Make_PDB_traj_from_images_v2(int Num_images, int Num_atoms, char *Filename_PDB_traj)
{
  int im, config_index;
  FR_DAT fr;
  char Filename_input[1000];
  FILE *out;
  int Num_particles;
  PDB_ATOM *PDB_atoms;
  
  fr.natoms = Num_atoms;
  fr.x = (rvec*)malloc(Num_atoms * sizeof(rvec));

  Read_alloc_PDB_info(&Num_particles, &PDB_atoms);
  if(Num_particles != Num_atoms)
    {
      printf("Number of particles from 'PDB_INFO' does not match command line input.\n");
      exit(1);
    }

  out = fopen(Filename_PDB_traj, "w");
  for(im = 0 ; im < Num_images ; im++)
    {
      sprintf(Filename_input, "COORDS_IMAGE_%d", (im+1));
      Read_config_simple_3d(Num_atoms, &fr, Filename_input);
      config_index = im + 1;
      Write_config_to_a_PDB_trajectory_file(&fr, PDB_atoms, config_index, out);
    }
  fprintf(out, "END\n");
  fclose(out);

  free(fr.x);
  free(PDB_atoms);
}

