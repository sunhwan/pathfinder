/* Calculate force and potential for a one state elastic network model.*/

/* Last modified on May 20, 2013*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>


#include "force_pot_simple_enm.h"

/* Calculate simple distance, square of distance and the distance vector 
 * between two points in d-dimension. *point_1 and *point_2 are 
 * d-diemnsional arrays (or d-dimensional vectors). This is meant for non-periodic systems
 * The distance_vector = point_2 - point_1*/
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



/* Calculate pair distances for both the structures and store them in separate arrays*/
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




/* calcualte force and potential*/
void getforces_enm(FR_DAT *fr, double **pair_distances_structure, double cutoff, double force_constant)
{
  int i,j;
  double distance, distance_squared;
  double distance_vector[3];
  double unit_vector[3];
  double temp_force[3];  /*this is used for various stages of force calculation*/
  double scalar_force, potential, diff;

  /* Initialize forces to zeros*/
  for(i = 0 ; i < fr->natoms ; i++)
    {
      fr->f[i][0] = 0.0;
      fr->f[i][1] = 0.0;
      fr->f[i][2] = 0.0;
    }
  
  /* Zero total energy*/
  fr->U = 0.0;
  
  /* Double loop over all possible pairs*/
  for(i = 0 ; i < (fr->natoms - 1) ; i++)
    {
      for(j = (i+1) ; j < fr->natoms ; j++)
	{
	  /* Calculate interaction if the pair is within interacting dstance in the native structure*/
          
	  if(pair_distances_structure[i][j] < cutoff)
	    {
	      /* Calculate distance and unit vector*/
	      Calc_simple_distance(3, fr->x[i], fr->x[j], &distance, &distance_squared, distance_vector);
              
	      unit_vector[0] = distance_vector[0] / distance;
	      unit_vector[1] = distance_vector[1] / distance;
	      unit_vector[2] = distance_vector[2] / distance;
	      diff = distance - pair_distances_structure[i][j];
	      potential = 0.5 * force_constant * pow(diff, 2);
	      fr->U = fr->U + potential;
	      scalar_force = -force_constant * diff;
	      temp_force[0] = scalar_force * unit_vector[0];
	      temp_force[1] = scalar_force * unit_vector[1];
	      temp_force[2] = scalar_force * unit_vector[2];
	      /* Sum force on i*/
	      fr->f[i][0] = fr->f[i][0] + (-temp_force[0]);
	      fr->f[i][1] = fr->f[i][1] + (-temp_force[1]);
	      fr->f[i][2] = fr->f[i][2] + (-temp_force[2]);
	      /* Sum force on j*/
	      fr->f[j][0] = fr->f[j][0] + temp_force[0];
	      fr->f[j][1] = fr->f[j][1] + temp_force[1];
	      fr->f[j][2] = fr->f[j][2] + temp_force[2];
	    }
	}
    }
}




/* Calculate force and potential. Return potential*/
double getforces_enm_return_pot(FR_DAT *fr, double **pair_distances_structure, double cutoff, double force_constant)
{
  int i,j;
  double distance, distance_squared;
  double distance_vector[3];
  double unit_vector[3];
  double temp_force[3];  /*this is used for various stages of force calculation*/
  double scalar_force, potential, diff;
  
  /* Initialize forces to zeros*/
  for(i = 0 ; i < fr->natoms ; i++)
    {
      fr->f[i][0] = 0.0;
      fr->f[i][1] = 0.0;
      fr->f[i][2] = 0.0;
    }
  
  /* Zero total energy*/
  fr->U = 0.0;
  
  /* Double loop over all possible pairs*/
  for(i = 0 ; i < (fr->natoms - 1) ; i++)
    {
      for(j = (i+1) ; j < fr->natoms ; j++)
	{
	  /* Calculate interaction if the pair is within interacting dstance in the native structure*/
	  
	  if(pair_distances_structure[i][j] < cutoff)
	    {
	      /* Calculate distance and unit vector*/
	      Calc_simple_distance(3, fr->x[i], fr->x[j], &distance, &distance_squared, distance_vector);
	      
	      unit_vector[0] = distance_vector[0] / distance;
	      unit_vector[1] = distance_vector[1] / distance;
	      unit_vector[2] = distance_vector[2] / distance;
	      diff = distance - pair_distances_structure[i][j];
	      potential = 0.5 * force_constant * pow(diff, 2);
	      fr->U = fr->U + potential;
	      scalar_force = -force_constant * diff;
	      temp_force[0] = scalar_force * unit_vector[0];
	      temp_force[1] = scalar_force * unit_vector[1];
	      temp_force[2] = scalar_force * unit_vector[2];
	      /* Sum force on i*/
	      fr->f[i][0] = fr->f[i][0] + (-temp_force[0]);
	      fr->f[i][1] = fr->f[i][1] + (-temp_force[1]);
	      fr->f[i][2] = fr->f[i][2] + (-temp_force[2]);
	      /* Sum force on j*/
	      fr->f[j][0] = fr->f[j][0] + temp_force[0];
	      fr->f[j][1] = fr->f[j][1] + temp_force[1];
	      fr->f[j][2] = fr->f[j][2] + temp_force[2];
	    }
	}
    }

  return fr->U;
}



/* calcualte force and potential. adds offset value to the potential*/
void getforces_enm_offset(FR_DAT *fr, double **pair_distances_structure, double cutoff, double force_constant, double energy_offset)
{
  int i,j;
  double distance, distance_squared;
  double distance_vector[3];
  double unit_vector[3];
  double temp_force[3];  /*this is used for various stages of force calculation*/
  double scalar_force, potential, diff;

  /* Initialize forces to zeros*/
  for(i = 0 ; i < fr->natoms ; i++)
    {
      fr->f[i][0] = 0.0;
      fr->f[i][1] = 0.0;
      fr->f[i][2] = 0.0;
    }
  
  /* Zero total energy*/
  fr->U = 0.0;
  
  /* Double loop over all possible pairs*/
  for(i = 0 ; i < (fr->natoms - 1) ; i++)
    {
      for(j = (i+1) ; j < fr->natoms ; j++)
	{
	  /* Calculate interaction if the pair is within interacting dstance in the native structure*/
          
	  if(pair_distances_structure[i][j] < cutoff)
	    {
	      /* Calculate distance and unit vector*/
	      Calc_simple_distance(3, fr->x[i], fr->x[j], &distance, &distance_squared, distance_vector);
              
	      unit_vector[0] = distance_vector[0] / distance;
	      unit_vector[1] = distance_vector[1] / distance;
	      unit_vector[2] = distance_vector[2] / distance;
	      diff = distance - pair_distances_structure[i][j];
	      potential = 0.5 * force_constant * pow(diff, 2);
	      fr->U = fr->U + potential;
	      scalar_force = -force_constant * diff;
	      temp_force[0] = scalar_force * unit_vector[0];
	      temp_force[1] = scalar_force * unit_vector[1];
	      temp_force[2] = scalar_force * unit_vector[2];
	      /* Sum force on i*/
	      fr->f[i][0] = fr->f[i][0] + (-temp_force[0]);
	      fr->f[i][1] = fr->f[i][1] + (-temp_force[1]);
	      fr->f[i][2] = fr->f[i][2] + (-temp_force[2]);
	      /* Sum force on j*/
	      fr->f[j][0] = fr->f[j][0] + temp_force[0];
	      fr->f[j][1] = fr->f[j][1] + temp_force[1];
	      fr->f[j][2] = fr->f[j][2] + temp_force[2];
	    }
	}
    }
  /* Add offset value to the energy*/
  fr->U = fr->U + energy_offset;
}




/* Calculate force and potential. Return potential. Add offset value to the potential*/
double getforces_enm_return_pot_offset(FR_DAT *fr, double **pair_distances_structure, double cutoff, double force_constant, double energy_offset)
{
  int i,j;
  double distance, distance_squared;
  double distance_vector[3];
  double unit_vector[3];
  double temp_force[3];  /*this is used for various stages of force calculation*/
  double scalar_force, potential, diff;
  
  /* Initialize forces to zeros*/
  for(i = 0 ; i < fr->natoms ; i++)
    {
      fr->f[i][0] = 0.0;
      fr->f[i][1] = 0.0;
      fr->f[i][2] = 0.0;
    }
  
  /* Zero total energy*/
  fr->U = 0.0;
  
  /* Double loop over all possible pairs*/
  for(i = 0 ; i < (fr->natoms - 1) ; i++)
    {
      for(j = (i+1) ; j < fr->natoms ; j++)
	{
	  /* Calculate interaction if the pair is within interacting dstance in the native structure*/
	  
	  if(pair_distances_structure[i][j] < cutoff)
	    {
	      /* Calculate distance and unit vector*/
	      Calc_simple_distance(3, fr->x[i], fr->x[j], &distance, &distance_squared, distance_vector);
	      
	      unit_vector[0] = distance_vector[0] / distance;
	      unit_vector[1] = distance_vector[1] / distance;
	      unit_vector[2] = distance_vector[2] / distance;
	      diff = distance - pair_distances_structure[i][j];
	      potential = 0.5 * force_constant * pow(diff, 2);
	      fr->U = fr->U + potential;
	      scalar_force = -force_constant * diff;
	      temp_force[0] = scalar_force * unit_vector[0];
	      temp_force[1] = scalar_force * unit_vector[1];
	      temp_force[2] = scalar_force * unit_vector[2];
	      /* Sum force on i*/
	      fr->f[i][0] = fr->f[i][0] + (-temp_force[0]);
	      fr->f[i][1] = fr->f[i][1] + (-temp_force[1]);
	      fr->f[i][2] = fr->f[i][2] + (-temp_force[2]);
	      /* Sum force on j*/
	      fr->f[j][0] = fr->f[j][0] + temp_force[0];
	      fr->f[j][1] = fr->f[j][1] + temp_force[1];
	      fr->f[j][2] = fr->f[j][2] + temp_force[2];
	    }
	}
    }

  /* Add offset value to the energy*/
  fr->U = fr->U + energy_offset;
  
  return fr->U;
}
