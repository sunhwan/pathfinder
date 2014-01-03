/* Align all the strucutres in the pathway with a reference strucutre. Then
 * make a PDB trajectory file out of aligned strucutres. Input files for images
 * are overwritten. Alignment is optional, the option is specified by a 
 * command line argument.*/

/* Uses new alignment code that can handle alignment based on a subset of atoms.*/
/* Uses 'PDB_INFO' as an input file. This is more general version of the previous
 * implementation.*/

/* Residue numbers, as present in the PDB_INFO file, are read as inputs.
   Atom indices are inferred from there.*/

/* Last updated on June 12, 2013*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "include_all.h"

/* 
 * Input file for system infomation to be used for writing PDB trajectory file: 'PDB_INFO'
 * Input files for coordinates: 'COORDS_IMAGE_1', 'COORDS_IMAGE_2', etc.
 *
 * Command line inputs: 1) Number of atoms.
 *                      2) Number of images.
 * Number of files = Number of images.
 *
 * Input file for optional alignment: 'REFERENCE_FOR_ALIGNMENT' and 'RESIDUES_FOR_ALIGNMENT'
 * The format of the 'RESIDUES_FOR_ALIGNMENT' file
 * -----------------------------------------------
 *   First chain id then either two integers per line separated by white space or one integer.
 *   Two integers indicate all the integers in the range will be included
 *   including the two. These integers are 1-based residue numbers present
 *   in the PDB_INFO file. 
 * -----------------------------------------------
 * Example
 * -------------------
 *  A 12
 *  B 13 56
 *  B 78
 *  C 113 200
 * -------------------
 *
 * The output file is called 'pathway.pdb'
 * */

/* Usage: ./makepathway_v2 -n Num_atoms -i Num_images -a Flag_indicating_alignment
          Flag_indicating_alignment = 1 if alignemnt is needed in that case files 'REFERENCE_FOR_ALIGNMENT' 
                                        and 'RESIDUES_FOR_ALIGNMENT' need to be present in the same folder.
                                    = 0 if alignment is not needed.
 * */


/* Returns the 1-based atom index for given 1-based residue number. Parses the 
 * 'PDB_INFO' file.*/
int Find_one_based_atom_index_from_one_based_residue_id(int Num_particles, PDB_ATOM *PDB_atoms, char chain_id, int residue_number_one_based)
{
  int i;
  for(i = 0 ; i < Num_particles ; i++)
    {
      if((PDB_atoms[i].chain_id == chain_id) && (PDB_atoms[i].residue_sequence_number == residue_number_one_based))
        {
          return (i + 1);
        }
    }
  /* If does not find it, returns -1*/
  return -1;
}

/* Reads the 'RESIDUES_FOR_ALIGNMENT' file.
   After reading it finds out the 1-based atom indices and writes 'SUBSET_FOR_ALIGNMENT' file. This file is read
   and used for alignment.
   Format for 'SUBSET_FOR_ALIGNMENT'
   ---------------------------------
    Number of atoms in the subset
    Atom indices one per line. Indices are 1-based meaning they start from 1.
   ---------------------------------
 */
void Read_residues_for_alignment(int Num_atoms, PDB_ATOM *PDB_atoms)
{
  FILE *in, *out;
  char line[10000];
  char chain_id;
  int resid_1, resid_2;
  int line_number;
  int *list_of_residues;     // List of residues to be included in the alignment. Read from 'RESIDUES_FOR_ALIGNMENT'
  int Num_residues;          // NUmber of residues to be included in the alignment.
  char *list_of_chain_ids;
  int r;
  int *list_of_one_based_indices;

  /* Maximum length of this array can be total number of atoms present in the system*/
  list_of_residues = (int*)malloc(Num_atoms * sizeof(int));
  list_of_chain_ids = (char*)malloc(Num_atoms * sizeof(char));

  in = fopen("RESIDUES_FOR_ALIGNMENT", "r");
  if(in == NULL)
    {
      printf("File 'RESIDUES_FOR_ALIGNMENT' not found\n");
      exit(1);
    }
  line_number = 1;
  Num_residues = 0;
  while(fgets(line, 10000, in) != NULL)
    {
      if(sscanf(line, "%c %d %d", &chain_id, &resid_1, &resid_2) == 3)
        {
          if(resid_2 < resid_1)
            {
              printf("First number should be smaller than the second one. Error in line %d of file 'RESIDUES_FOR_ALIGNMENT'.\n", line_number);
              exit(1);
            }
          for(r = resid_1 ; r <= resid_2 ; r++)
            {
              list_of_residues[Num_residues] = r;
              list_of_chain_ids[Num_residues] = chain_id;
              Num_residues = Num_residues + 1;
            }
        }
      else if(sscanf(line, "%c %d", &chain_id, &resid_1) == 2)
        {
          list_of_residues[Num_residues] = resid_1;
          list_of_chain_ids[Num_residues] = chain_id;
          Num_residues = Num_residues + 1;
        }
      else
        {
          printf("Error in line %d of file 'RESIDUES_FOR_ALIGNMENT'.\n", line_number);
        }
      line_number = line_number + 1;
    }
  fclose(in);

  list_of_one_based_indices = (int*)malloc(Num_residues * sizeof(int));
  /* Dtermine 1-based particle indices from the 1-based residue indices*/
  for(r = 0 ; r < Num_residues ; r++)
    {
      list_of_one_based_indices[r] = Find_one_based_atom_index_from_one_based_residue_id(Num_atoms, PDB_atoms, list_of_chain_ids[r], list_of_residues[r]);
      if(list_of_one_based_indices[r] == -1)
        {
          printf("Residue id %d with chain id %c not found in the PDB_INFO file.\n", list_of_residues[r], list_of_chain_ids[r]);
          exit(1);
        }
    }

  /* Write 'SUBSET_FOR_ALIGNMENT' file*/ 
  out = fopen("SUBSET_FOR_ALIGNMENT", "w");
  fprintf(out, "%d\n", Num_residues);
  for(r = 0 ; r < Num_residues ; r++)
    {
      fprintf(out, "%d\n", list_of_one_based_indices[r]);
    }
  fclose(out);  

  free(list_of_residues);
  free(list_of_one_based_indices);
}


int main(int argc, char *argv[])
{
  char *error_message = "Usage: ./makepathway_v2 -n Num_atoms -i Num_images -a Flag_indicating_alignment\nFlag_indicating_alignment = 1 if alignemnt is needed, in that case files 'REFERENCE_FOR_ALIGNMENT'\n                              and 'RESIDUES_FOR_ALIGNMENT' need to be present in the same folder.\n                          = 0 if alignment is not needed.\n";
  if(argc == 1)
    {
      printf("For help type: ./makepathway_v2 -h\n");
      exit(1);
    }
  if(strcmp(argv[1],"-h")==0)
    {
      //printf("The usage is: ./makepathway -n Num_atoms -i Num_images\n");
      printf("%s", error_message);
      exit(1);
    }
  if(argc != 7)
    {
      printf("Error!\n");
      //printf("The usage is: ./makepathway -n Num_atoms -i Num_images\n");
      printf("%s", error_message);
      exit(1);
    }
  if(strcmp(argv[1],"-n")!=0 || strcmp(argv[3],"-i")!=0 || strcmp(argv[5],"-a")!=0)
    {
      printf("Error!\n");
      //printf("The usage is: ./makepathway -n Num_atoms -i Num_images\n");
      printf("%s", error_message);
      exit(1);
    }

  int Num_images;
  int Num_atoms; 
  char Filename_PDB_traj[1000];
  int Num_atoms_in_subset;
  int *Indices_atoms_in_subset;
  FR_DAT fr_ref, fr;
  char Filename[1000];
  int im;
  int align_flag;
  PDB_ATOM *PDB_atoms;
  int Num_particles_from_PDB_INFO;

  Num_atoms = atoi(argv[2]);
  Num_images = atoi(argv[4]);
  align_flag = atoi(argv[6]);
  if(align_flag != 0 && align_flag != 1)
    {
      printf("Input error.\n -a should be followed by either 1 or 0.\n");
      exit(1);
    }

  fr_ref.natoms = Num_atoms;
  fr_ref.x = (rvec*)malloc(Num_atoms * sizeof(rvec));
  fr.natoms = Num_atoms;
  fr.x = (rvec*)malloc(Num_atoms * sizeof(rvec));

  /* Check if alignment is requested or not.*/
  if(align_flag == 1)
    {
      /* Read the reference structure*/
      sprintf(Filename, "REFERENCE_FOR_ALIGNMENT");
      Read_config_simple_3d(Num_atoms, &fr_ref, Filename);

      /* Read PDB_INFO*/
      Read_alloc_PDB_info(&Num_particles_from_PDB_INFO, &PDB_atoms);
      if(Num_particles_from_PDB_INFO != Num_atoms)
        {
          printf("Inconsistent input for number of particles\n");
          exit(1);
        }

      /* Read residues for alignment*/
      Read_residues_for_alignment(Num_atoms, PDB_atoms);
      
      /* Read information for alignment*/
      Read_alloc_subset_for_alignment(Num_atoms, &Num_atoms_in_subset, &Indices_atoms_in_subset);
      
      /* Align all strucutres in the pathway*/
      for(im = 0 ; im < Num_images ; im++)
      	{
      	  /* Read image*/
      	  sprintf(Filename, "COORDS_IMAGE_%d", (im+1));
      	  Read_config_simple_3d(Num_atoms, &fr, Filename);
	  /* Align image*/
      	  Align_two_strucutres_using_subset_lei(&fr_ref, Num_atoms_in_subset, Indices_atoms_in_subset, &fr);
	  /* Write image, overwrite the input file*/
      	  Write_config_simple(&fr, Filename);
      	}
    }
  
  /* Make PDB trajectory file from aligned images*/
  sprintf(Filename_PDB_traj, "pathway.pdb");
  Make_PDB_traj_from_images_v2(Num_images, Num_atoms, Filename_PDB_traj);
  
  return 0;
}


