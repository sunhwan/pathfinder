/* Codies for writing PDB files with one cofromation, multiple conformation etc.
 * These functions are more general than what we had before, the present codes can handle
 * files with multiple chains, missing residues etc.
 * Different type of input information is employed.
 * See below for details.
 * */

/* last updated July 09, 2012*/




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "write_PDB_files.h"


/* This function reads some information about the system from the 'PDB_INFO' file.
 * The content of this file is the minimum amount of information needed to specify the system.
 * There are are other less important infomration in a PDB file, those are not stored in 
 * the 'PDB_INFO' file. Given a set of coordinates and the 'PDB_INFO' file a complete PDB 
 * file can be generated after populating some less important fields with generic values.
 * The format of 'PDB_INFO' file:
 * ---------------------------------
 *  Number of atoms (%d)
 *  Each line will have the follwoing fields separated by white space
 *    (Atom_serial_number) (Atom_name) (Three_letter_residue_name) (Chain_identifier) (Residue_serial_number) (Element_symbol)
 * ---------------------------------
 * These fields are extracted from the PDB file, if entered from scrath then they have to 
 * be consistent with the official PDB format.
 **/
void Read_alloc_PDB_info(int *Num_particles, PDB_ATOM **PDB_atoms)
{
  FILE *in;
  int i, line_number;

  in = fopen("PDB_INFO", "r");
  if(in == NULL)
    {
      printf("File 'PDB_INFO' not found\n");
      exit(1);
    }
  line_number = 1;
  if(fscanf(in, "%d ", Num_particles) != 1)
    {
      printf("Error in 'PDB_INFO' file. Line %d.\n", line_number);
      exit(1);
    }
  line_number++;
  /* Allocate memory*/
  *PDB_atoms = (PDB_ATOM*)malloc((*Num_particles) * sizeof(PDB_ATOM));
  if(*PDB_atoms == NULL)
    {
      printf("Memory allocation error in 'write_PDB_files.c'. 1\n");
      exit(1);
    }
  for(i = 0 ; i < *Num_particles ; i++)
    {
      if(fscanf(in, "%d %s %s %c %d %s ", 
		&(*PDB_atoms)[i].atom_serial_number,
		(*PDB_atoms)[i].atom_name,
		(*PDB_atoms)[i].residue_name,
		&(*PDB_atoms)[i].chain_id,
		&(*PDB_atoms)[i].residue_sequence_number,
                (*PDB_atoms)[i].element_symbol) != 6)
	{
	  printf("Error in 'PDB_INFO' file. Line %d.\n", line_number);
	  exit(1);
	}
      line_number++;
    }
  fclose(in);
}


/* Writes a frame or configuration to a PDB file. Input: coordinates and essential info. Some
 * fields are populated with generic values.*/
void Write_one_config_to_PDB_file(FR_DAT *fr, PDB_ATOM *PDB_atoms, FILE *out)
{
  int i;
 
  for(i = 0 ; i < fr->natoms ; i++)
    {
      /* Populate empty fields, some with generic values.*/
      PDB_atoms[i].atom_alt_loc = ' ';
      PDB_atoms[i].insert_code = ' ';
      /*need to typecast positions values from double to float*/
      PDB_atoms[i].x = (float)fr->x[i][0];
      PDB_atoms[i].y = (float)fr->x[i][1];
      PDB_atoms[i].z = (float)fr->x[i][2];
      PDB_atoms[i].occupancy = 1.00;
      PDB_atoms[i].temperature_factor = 1.00;
      sprintf(PDB_atoms[i].charge_on_atom, " ");

      /* Write to file*/
      Write_one_PDB_ATOM(out, &PDB_atoms[i]); 
    }



}


/* Writes a PDB file*/
void Write_PDB_FILE(FR_DAT *fr, PDB_ATOM *PDB_atoms, char *Filename_PDB)
{
  FILE *out;
  
  out = fopen(Filename_PDB, "w");
  Write_one_config_to_PDB_file(fr, PDB_atoms, out);
  fprintf(out, "END\n");
  fclose(out);
}




/* This function writes a configuration to a PDB trajectory file.
 * 'config_index' gives the index of the configuration written. It starts from 1 and takes values upto total Number of
 * trajectories in the file.
 *
 * The format of the file. Each config is written in the following way.
 * MODEL        config_index
 * Info about config in PDB format
 * ENDMDL
 *
 * */
void Write_config_to_a_PDB_trajectory_file(FR_DAT *fr, PDB_ATOM *PDB_atoms, int config_index, FILE *out)
{
  fprintf(out, "MODEL        %d\n", config_index);
  Write_one_config_to_PDB_file(fr, PDB_atoms, out); 
  fprintf(out, "ENDMDL\n");
}



