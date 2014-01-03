/* Find pairs that come closer during the transition but are farther apart 
   in the end states
 */

/* Reads PDB_INFO file
   Outputs the residue numbers of contact forming pairs instead the atom numbers.
   Distance profiles files are also indexed by residue numbers instead of ato mnumbers.*/

/* Last modified on August 23, 2013*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "include_all.h"


/* Input filename: 'INPUT_FIND_CLOSE_PAIRS'
 * Format:
 * ---------------------------
 *  NUM_PARTICLES (Number of particles in the system. %d)
 *  NUM_IMAGES (Number of images in the pathway. %d)
 *  CUT_OFF_PATHWAY (Cut-off distance for the pathway %lf)
 *  OTHER_DISTANCE (Other distance for the end states %lf)
 * ---------------------------
 * */


void Read_inputs(int *Num_particles, int *Num_images, double *cut_off_pathway, double *separation_end_states)
{
  FILE *in;
  char key[1000];

  in = fopen("INPUT_FIND_CLOSE_PAIRS", "r");
  if(in == NULL)
    {
      printf("File 'INPUT_FIND_CLOSE_PAIRS' not found\n");
      exit(1);
    }
  fscanf(in, "%s  %d ", key, Num_particles);
  fscanf(in, "%s %d ", key, Num_images);
  fscanf(in, "%s %lf ", key, cut_off_pathway);
  fscanf(in, "%s %lf ", key, separation_end_states);
  fclose(in);
}


/* Find minimum of two numbers.*/
double Find_minimum(double number_1, double number_2)
{
  double min_value;
  if(number_1 <= number_2)
    {
      min_value = number_1;
    }
  else
    {
      min_value = number_2;
    }

  return min_value;
}


void Calc_simple_distance(int dimension, double *point_1, double *point_2, double *distance, double *distance_squared, double *distance_vector)
{
  int i;
  for(i=0;i<dimension;i++)
    {
      distance_vector[i] = point_2[i] - point_1[i];
    }
  *distance_squared = 0.0;
  for(i=0;i<dimension;i++)
    {
      *distance_squared = *distance_squared + (distance_vector[i] * distance_vector[i]);
    }
  *distance = sqrt(*distance_squared);
}



void Calc_pair_distances_both_structures(int Num_particles, FR_DAT fr_structure_1, FR_DAT fr_structure_2, double **pair_distances_structure_1, double **pair_distances_structure_2)
{
  int i, j;
  int dimension = 3;
  double distance;
  double distance_squared;
  double distance_vector[3];

  /* Set diagonal elements to -1.0*/
  for(i = 0 ; i < Num_particles ; i++)
    {
      pair_distances_structure_1[i][i] = -1.0;
      pair_distances_structure_2[i][i] = -1.0;
    }

  for(i = 0 ; i < (Num_particles - 1) ; i++)
    {
      for(j = (i + 1) ; j < Num_particles ; j++)
        {
          /* Structure 1*/
          Calc_simple_distance(dimension, fr_structure_1.x[i], fr_structure_1.x[j], &distance, &distance_squared, distance_vector);
          pair_distances_structure_1[i][j] = distance;
          pair_distances_structure_1[j][i] = distance;

          /* Structure 2*/
          Calc_simple_distance(dimension, fr_structure_2.x[i], fr_structure_2.x[j], &distance, &distance_squared, distance_vector);
          pair_distances_structure_2[i][j] = distance;
          pair_distances_structure_2[j][i] = distance;
        }
    }

}



/* Calculate pair distances for one end strucutre*/
void Calc_pair_distances_one_structure(int Num_particles, FR_DAT fr_structure, double **pair_distances_structure)
{
  int i, j;
  int dimension = 3;
  double distance;
  double distance_squared;
  double distance_vector[3];

  /* Set diagonal elements to -1.0*/
  for(i = 0 ; i < Num_particles ; i++)
    {
      pair_distances_structure[i][i] = -1.0;
    }

  for(i = 0 ; i < (Num_particles - 1) ; i++)
    {
      for(j = (i + 1) ; j < Num_particles ; j++)
        {
          /* Structure*/
          Calc_simple_distance(dimension, fr_structure.x[i], fr_structure.x[j], &distance, &distance_squared, distance_vector);
          pair_distances_structure[i][j] = distance;
          pair_distances_structure[j][i] = distance;
        }
    }
}


/* Finds 1-based residue number and chain id given the 1-based atom index. Parses the 
   content of the 'PDB_INFO' file */
void Find_chain_id_res_num_from_1_based_atom_index(int Num_particles, PDB_ATOM *PDB_atoms, int atom_index_1_based, char *chain_id, int *residue_number_one_based)
{
  int i;

  i = atom_index_1_based - 1;
  *chain_id = PDB_atoms[i].chain_id;
  *residue_number_one_based = PDB_atoms[i].residue_sequence_number;
}


int main()
{
  int Num_particles;
  int Num_images;
  double cut_off_pathway;
  double separation_end_states;
  FR_DAT fr_ref_1, fr_ref_2, fr;
  double **pair_distances_structure_1, **pair_distances_structure_2, **pair_distances_current;
  double **min_separation_end;
  int **final_pairs;
  char Filename[1000];
  int i, j, im;
  double min_dis;
  double distance, distance_squared;
  double distance_vector[3];
  int dimension;
  PDB_ATOM *PDB_atoms;
  int Num_particles_from_PDB_INFO;
  int i_1, j_1;
  int residue_number_i, residue_number_j;
  char chain_id_i, chain_id_j;

  /****************************************************************************************/
  
  dimension = 3;

  Read_inputs(&Num_particles, &Num_images, &cut_off_pathway, &separation_end_states);

  Read_alloc_PDB_info(&Num_particles_from_PDB_INFO, &PDB_atoms);
  if(Num_particles_from_PDB_INFO != Num_particles)
    {
      printf("Inconsistent input for number of particles\n");
      exit(1);
    }

  pair_distances_structure_1 = (double**)malloc(Num_particles * sizeof(double*));
  pair_distances_structure_2 = (double**)malloc(Num_particles * sizeof(double*));
  pair_distances_current = (double**)malloc(Num_particles * sizeof(double*));
  min_separation_end = (double**)malloc(Num_particles * sizeof(double*));
  final_pairs = (int**)malloc(Num_particles * sizeof(int*));
  for(i = 0 ; i < Num_particles ; i++)
    {
      pair_distances_structure_1[i] = (double*)malloc(Num_particles * sizeof(double));
      pair_distances_structure_2[i] = (double*)malloc(Num_particles * sizeof(double));
      pair_distances_current[i] = (double*)malloc(Num_particles * sizeof(double));
      min_separation_end[i] = (double*)malloc(Num_particles * sizeof(double));
      final_pairs[i] = (int*)malloc(Num_particles * sizeof(int));
    }

  for(i = 0 ; i < Num_particles ; i++)
    {
      for(j = 0 ; j < Num_particles ; j++)
        {
          final_pairs[i][j] = 0;
        }
    }

  fr_ref_1.natoms = Num_particles;
  fr_ref_1.x = (rvec*)malloc(Num_particles * sizeof(rvec));
  fr_ref_2.natoms = Num_particles;
  fr_ref_2.x = (rvec*)malloc(Num_particles * sizeof(rvec));

  fr.natoms = Num_particles;
  fr.x = (rvec*)malloc(Num_particles * sizeof(rvec));


  sprintf(Filename, "INPUT_STRUCTURE_1");
  Read_config_simple_3d(Num_particles, &fr_ref_1, Filename);
  sprintf(Filename, "INPUT_STRUCTURE_2");
  Read_config_simple_3d(Num_particles, &fr_ref_2, Filename);
 
  Calc_pair_distances_both_structures(Num_particles, fr_ref_1, fr_ref_2, pair_distances_structure_1, pair_distances_structure_2);

  for(i = 0 ; i < (Num_particles - 1) ; i++)
    {
      for(j = i+1 ; j < Num_particles ; j++)
        {
          min_dis = Find_minimum(pair_distances_structure_1[i][j], pair_distances_structure_2[i][j]);
          min_separation_end[i][j] = min_dis;
          min_separation_end[j][i] = min_dis; 
        } 
    }


  for(im = 0 ; im < Num_images ; im++)
    {
      sprintf(Filename, "COORDS_IMAGE_%d", (im+1));
      Read_config_simple_3d(Num_particles, &fr, Filename);
      Calc_pair_distances_one_structure(Num_particles, fr, pair_distances_current);
      for(i = 0 ; i < (Num_particles - 1) ; i++)
        {
          for(j = i+1 ; j < Num_particles ; j++)
            {
              if((pair_distances_current[i][j] < cut_off_pathway) && (min_separation_end[i][j] > separation_end_states))
                {
                  final_pairs[i][j] = 1;
                  final_pairs[j][i] = 1;
                }
            }
        }
    }
  

  FILE *out, *out_1;
  char Filename_1[1000], Filename_2[1000];
  sprintf(Filename, "close_contacts_%1.1f_%2.1f", cut_off_pathway, separation_end_states);
  out = fopen(Filename, "w");
  for(i = 0 ; i < (Num_particles - 1) ; i++)
    {
      for(j = i+1 ; j < Num_particles ; j++)
        {
          if(final_pairs[i][j] == 1)
            {
              i_1 = i+1;
              j_1 = j+1;
              Find_chain_id_res_num_from_1_based_atom_index(Num_particles, PDB_atoms, i_1, &chain_id_i, &residue_number_i); 
              Find_chain_id_res_num_from_1_based_atom_index(Num_particles, PDB_atoms, j_1, &chain_id_j, &residue_number_j);
              fprintf(out, "%c %d      %c %d\n", chain_id_i, residue_number_i, chain_id_j, residue_number_j);
              sprintf(Filename_1, "distance_%c%d_%c%d_%1.1f_%2.1f", chain_id_i, residue_number_i, chain_id_j, residue_number_j, cut_off_pathway, separation_end_states);
              out_1 = fopen(Filename_1, "w");
              for(im = 0 ; im < Num_images ; im++)
                {
                  sprintf(Filename_2, "COORDS_IMAGE_%d", (im+1));
                  Read_config_simple_3d(Num_particles, &fr, Filename_2);
                  Calc_simple_distance(dimension, fr.x[i], fr.x[j], &distance, &distance_squared, distance_vector);
                  fprintf(out_1, "%d   %lf\n", (im+1), distance);
                }
              fclose(out_1);
            }
        }
    }
  fclose(out);

  return 0;
}




