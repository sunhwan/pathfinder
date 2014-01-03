/* This file contains common data structures for performing 
 * molecualr dynamics (MD) or Langevin dynamics (LD) 
 * on a multiple component d-dimensional molecular system.
   No electrostatic interactions, no intramolecular constraints.
*/

/* To change the dimensionality of the system change the value of DIMENSION.*/

/* These structures can also be used for simple systems like Muller potential.*/

/* For Langevin dynamics we will use the integrator proposed by Vanden-Eijnden and 
 * Ciccotti Chem. Phys. Lett. vol 429, pp 310-316, 2006
 * 
 * But these structures can presumably be used with other integrators also.*/

/* This structures are analogues of GROMACS strucutres and some structures 
   used by Lanyuan Lu for his force matching program. But some modifications have been used later.*/

/* We change the datatype 'matrix' to 'my_matrix'. Because Numerical Recipes has its own matrix
 * data type and there was a conflict in namespace while compilation. Pay attention to this 
 * when using this header file for other purposes.*/

/* NAMD units are used. See below*/

/* Units:
 * Internal
 * --------
 *  Energy = kcal mol^{-1}
 *  Mass = amu
 *  Distance = \AA (angstrom)
 *  Time = 48.88821291 fs (femto second)
 * 
 * External
 * --------
 *  Energy = kcal mol^{-1}
 *  Mass = amu
 *  Distance = \AA (angstrom)
 *  Time = fs (femto second)
 * */

/* External units are used for all input output for positions and times.
 * I.e. any quantity that has dimension of length should be in angstrom in both input and output
 * operations. Similarly for times femtosecond is used for both input and output.
 * */

/* Unit conversion for qunatities with dimension of time (e.g. timestep) is handled internally.*/

/* IF NEW UNITS ARE INTRODUCED THEN VALUE OF BOLTZMANN CONSTANT SHOULD BE MODIFIED AS WELL AS
 * UNIT OF TIME SHOULD BE WORKED OUT PROPERLY. THIS MEANS MACROS THAT CONVERT BETWEEN TIME UNITS
 * (FOR THE PRESENT VERSION THESE ARE 'INTERNAL_TIME_TO_FS' AND 'FS_TO_INTERNAL_TIME') WILL HAVE
 * DIFFERENT VALUES.*/

/* Written by Avisek Das. */

/* Last modified on December 05, 2011 */

#ifndef _md_ld_common_d_dimension_namd_h
#define _md_ld_common_d_dimension_namd_h


#define VERYSMALL    1.0e-14
#define VERYSMALL_F  1.0e-6          //small number for single precision
#define DIMENSION    3              /* Spatial diemnsion.Even though part of 
                                        this code is flexible to change of 
                                        dimension, most of it is not. It is 
                                        strictly applicable for three 
                                        dimensional systems. */
#define R2D          57.2957795      //radian to degree conversion factor
#define D2R          0.0174532925    //degree to radian conversion factor
/*when energy is in kJ mol^-1, temperature K. These are GROMACS units  */
//#define BOLTZMANN_CONSTANT   0.00831451
//#define BOLTZMANN_CONSTANT 1.0
/* when energy is in kcal mol^-1, temperature K. These are NAMD units */
#define BOLTZMANN_CONSTANT 0.0019872041 // in kcal mol^{-1} K^{-1}

#define INTERNAL_TIME_TO_FS 48.88821291   //conversion factor for internal time unit to femtosecond
#define FS_TO_INTERNAL_TIME 0.02045482828 //conversion factor from femtosecond to internal time unit
#define INVERSE_FS_TO_INVERSE_INTERNAL_TIME 48.88821291
#define INVERSE_INTERNAL_TIME_TO_INVERSE_FS 0.02045482828 


/* Data types. This is analogous to GROMACS data types.*/
typedef double my_matrix[DIMENSION][DIMENSION];
typedef double rvec[DIMENSION];


/*Definitions of data structures used*/



/* This is the data structure for storing the particle coordinates, forces and
   velocities for one frame (or configuration) of the trajectory. This is the main 
   data structure seen by all funtions in the MD code whenever the functions need 
   coordinate, force or velocity information.*/
typedef struct
{
  int natoms;         //total number of atoms
  my_matrix box;         //box dimension
  double box_dim2[DIMENSION];   //half box length

  rvec *x;          /*Coordinates of all the sites. 
	              If *fr is the pointer to the data structure 
	              for storing the information for one frame.
	              Then
	              x component of position of i th site = fr->x[i][0]
	              y component of position of i th site = fr->x[i][1]
	              z component of position of i th site = fr->x[i][2]
                      Straigtforward generalization for systems with other 
                      dimnesions.
	             */
  rvec *f;          /*Total forces on all the atoms.
	              Forces on individual atoms can be accessed in the 
	              same manner as the coordinates.
	             */
  rvec *v;          /*Velocities of all the sites.*/

  double *m;                /* Masses of the atoms */

  double temperature;             /* Temperature of the system*/

  double t;                /* time */
  double E;                /* total energy of the system*/
  double K;                /* total kinetic energy of the system*/
  double U;                /* total potential energy of the system*/
  double U_r;              /* energy from the purely replusive part of the potential*/
  double U_en;             /* energy from the elastic net part*/
  double P;                /* pressure */
  
  
  int ncell;    //number of cells
  int *list;    //linked-lists
  int *head;    //head of linked-lists
  int m0;       //number of cells in x-direction
  int m1;       //number of cells in y-direction
  int m2;       //number of cells in z-direction
  int *map;     //neighbor cell map
  double cellx;   //cell dimension in x-direction
  double celly;   //cell dimension in y-direction
  double cellz;   //cell dimension in z-direction
} FR_DAT;



/*Simulation details are stored in the structure 'simudetails'*/
struct simudetails
{
  long int RUNLENGTH;               /*Total run length in unit of timesteps. */
  double   timestep;                /*Timestep used in simulation*/
  int      ColliFreqInTimestep;     /*frequency of thermalising collision in 
				      the unit of time-steps*/
  int      NumFrames;               /*No. of frames for averaging g(r).*/
  double   MaxDistForRDF;           /*this is the maximum distance for which 
				      g(r) is calculated*/
  double   deltar;                  /*grid spacing for the distance grid used 
				      for g(r) calculation*/
  long int      ConfigSaveFreqTimeStep;  /*frequency in the unit of timestep at which 
				      configurations are saved to disk.*/
  long int      ProgressFreqTimeStep;    /*frequency of wrting progress info*/
  long int      WriteFreqTimestep;       /*frequency of writing other stuff to output files*/
  long int      WriteFreqRMSD;           /*frequency for writing RMSD to an output file.*/
  long int      WriteFreqPDBtraj;        /*frequency for writing configurations to a PDB trajectory file.*/
};


/* Inputs for Langevin dynamics.*/
struct langevin
{
  double *gamma;                   /*An array with dimension equal to the number of particles
                                    in the system that stores the friction coefficients
                                    of all the particles. The same value of friction coefficient
                                    is used for all the components of a single particle.
                                    I.e. for a three dimensional system 
                                    \gamma_{ix} = \gamma_{ix} = \gamma_{ix} = \gamma_i*/

  double *sigma_langevin;         /*An array of length equa to the number of particles in the system with
                                    the following values.

                                    \sigma_i = \sqrt(2 k_B T \gamma_i / m_i)

                                    where k_B is the Boltzmann constant, m_i is the mass of the i th 
                                    particle and \gamma_i is the friction coefficient for i th particle.*/
  
};


#endif
