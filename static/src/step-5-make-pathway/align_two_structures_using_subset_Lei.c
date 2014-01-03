/* Align two strucutres using the coordinates of a subset of atoms using Lei Huang's 
 * codes.*/

/* Last updated on July 09, 2012*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/*This is from Lei's code.*/
#define ZERO    1.0E-10

#include "align_two_structures_using_subset_Lei.h"



/* This is unchanged from Lei's code.*/
void jacobi (int n, int np, double a[5][5], double d[5], double v[5][5], double b[5], double z[5])
{
	int i,j,k;
	int ip,iq;
	int nrot,maxrot;
	double sm,tresh,s,c,t;
	double theta,tau,h,g,p;
	
	maxrot = 100;
	nrot = 0;
	
	for(ip=1; ip<=n; ip++)	{
		for(iq=1; iq<=n; iq++)	{
			v[ip][iq] = 0.0;
		}
		v[ip][ip] = 1.0;
	}
	for(ip=1; ip<=n; ip++)	{
		b[ip] = a[ip][ip];
		d[ip] = b[ip];
		z[ip] = 0.0;
	}
	
	
	for(i=1; i<=maxrot; i++)	{
		sm=0.0;
		for(ip=1; ip<=(n-1); ip++)	{
			for(iq=ip+1; iq<=n; iq++)	{
				sm = sm + fabs(a[ip][iq]);
			}
		}
		
		if(sm < ZERO)	{	//converge
			break;
		}
		if(i < 4)	{
			tresh = 0.2*sm / (n*n);
		}
		else	{
			tresh=0.0;
		}
		
		for(ip=1; ip<=(n-1); ip++)	{
			for(iq=ip+1; iq<=n; iq++)	{
				g = 100.0 * fabs(a[ip][iq]);
				if((i > 4) && (fabs(g) < ZERO))	{
					a[ip][iq] = 0.0;
				}
				else if(fabs(a[ip][iq]) > tresh)	{
					h = d[iq] - d[ip];
					if(fabs(g) < ZERO)	{
						t = a[ip][iq] / h;
					}
					else	{
						theta = 0.5*h / a[ip][iq];
						t = 1.0 / (fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0)  t = -t;
					}
					
					c = 1.0 / sqrt(1.0+t*t);
					s = t * c;
					tau = s / (1.0+c);
					h = t * a[ip][iq];
					z[ip] = z[ip] - h;
					z[iq] = z[iq] + h;
					d[ip] = d[ip] - h;
					d[iq] = d[iq] + h;
					a[ip][iq] = 0.0;
					
					for(j=1; j<=(ip-1); j++)	{
						g = a[j][ip];
						h = a[j][iq];
						a[j][ip] = g - s*(h+g*tau);
						a[j][iq] = h + s*(g-h*tau);
					}
					for(j=ip+1; j<=(iq-1); j++)	{
						g = a[ip][j];
						h = a[j][iq];
						a[ip][j] = g - s*(h+g*tau);
						a[j][iq] = h + s*(g-h*tau);
					}
					for(j=iq+1; j<=n; j++)	{
						g = a[ip][j];
						h = a[iq][j];
						a[ip][j] = g - s*(h+g*tau);
						a[iq][j] = h + s*(g-h*tau);
					}
					for(j=1; j<=n; j++)	{
						g = v[j][ip];
						h = v[j][iq];
						v[j][ip] = g - s*(h+g*tau);
						v[j][iq] = h + s*(g-h*tau);
					}
					nrot = nrot + 1;
				}
			}
		}
		for(ip=1; ip<=n; ip++)	{
            b[ip] = b[ip] + z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
		}
	}
	if(nrot >= maxrot)	{
		printf("JACOBI  --  Matrix Diagonalization not Converged.\n");
	}

	for(i=1; i<= (n-1); i++)	{
         k = i;
         p = d[i];
		 for(j=i+1; j<=n; j++)	{
			 if(d[j] < p)	{
               k = j;
               p = d[j];
			 }
		 }
		 if(k != i)	{
            d[k] = d[i];
            d[i] = p;
			for(j=1; j<=n; j++)	{
               p = v[j][i];
               v[j][i] = v[j][k];
               v[j][k] = p;
			}
		 }
	}
}




/* Calculate centers of mass. This is construcutred from Lei's code. All arrays have
 * indices starting from 1. The positions of centers of mass are stored.
 * In case of subset of atoms this is to be used with the subset not all the atoms*/
void Calc_center(double *xa, double *ya, double *za, double *xb, double *yb, double *zb, int AtomNum, double *midxa, double *midya, double *midza, double *midxb, double *midyb, double *midzb)
{
  int i;
  
  *midxa=*midya=*midza=*midxb=*midyb=*midzb=0.0; 
  
  for(i=1; i<=AtomNum; i++)       
    {
      *midxa += (xa[i]);
      *midya += (ya[i]);
      *midza += (za[i]);
      
      *midxb += (xb[i]);
      *midyb += (yb[i]);
      *midzb += (zb[i]);
    }
  *midxa/=AtomNum;
  *midya/=AtomNum;
  *midza/=AtomNum;
  *midxb/=AtomNum;
  *midyb/=AtomNum;
  *midzb/=AtomNum;

}

/* Translate to the respective center of mass. This is to be done with all the atoms*/
void Translate(double *xa, double *ya, double *za, double *xb, double *yb, double *zb, int AtomNum, double midxa, double midya, double midza, double midxb, double midyb, double midzb)
{
  int i;
  
  /* Translate*/
  for(i=1; i<=AtomNum; i++)       
    {
      xa[i]-=midxa;
      ya[i]-=midya;
      za[i]-=midza;
      xb[i]-=midxb;
      yb[i]-=midyb;
      zb[i]-=midzb;
    }
}


/* Calcualte the roation matrix. I ncase of subset this is to be calcualted using the subset only*/
void Calc_rotation_matrix(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n, double rot[4][4])
{
  double xxyx,xxyy,xxyz, xyyx,xyyy,xyyz, xzyx,xzyy,xzyz;
  double work1[5],work2[5], q[5],d[5], c[5][5], v[5][5];
  int ia;
  
  xxyx = 0.0;
  xxyy = 0.0;
  xxyz = 0.0;
  xyyx = 0.0;
  xyyy = 0.0;
  xyyz = 0.0;
  xzyx = 0.0;
  xzyy = 0.0;
  xzyz = 0.0;

  for(ia=1; ia<=n; ia++)  
    {
      xxyx = xxyx + x1[ia]*x2[ia];
      xxyy = xxyy + y1[ia]*x2[ia];
      xxyz = xxyz + z1[ia]*x2[ia];
      xyyx = xyyx + x1[ia]*y2[ia];
      xyyy = xyyy + y1[ia]*y2[ia];
      xyyz = xyyz + z1[ia]*y2[ia];
      xzyx = xzyx + x1[ia]*z2[ia];
      xzyy = xzyy + y1[ia]*z2[ia];
      xzyz = xzyz + z1[ia]*z2[ia];
    }

  c[1][1] = xxyx + xyyy + xzyz;
  c[1][2] = xzyy - xyyz;
  c[2][2] = xxyx - xyyy - xzyz;
  c[1][3] = xxyz - xzyx;
  c[2][3] = xxyy + xyyx;
  c[3][3] = xyyy - xzyz - xxyx;
  c[1][4] = xyyx - xxyy;
  c[2][4] = xzyx + xxyz;
  c[3][4] = xyyz + xzyy;
  c[4][4] = xzyz - xxyx - xyyy;
  
  jacobi(4,4,c,d,v,work1,work2);
  q[1] = v[1][4];
  q[2] = v[2][4];
  q[3] = v[3][4];
  q[4] = v[4][4];

  rot[1][1] = q[1]*q[1] + q[2]*q[2] - q[3]*q[3] - q[4]*q[4];
  rot[2][1] = 2.0 * (q[2] * q[3] - q[1] * q[4]);
  rot[3][1] = 2.0 * (q[2] * q[4] + q[1] * q[3]);
  rot[1][2] = 2.0 * (q[3] * q[2] + q[1] * q[4]);
  rot[2][2] = q[1]*q[1] - q[2]*q[2] + q[3]*q[3] - q[4]*q[4];
  rot[3][2] = 2.0 * (q[3] * q[4] - q[1] * q[2]);
  rot[1][3] = 2.0 * (q[4] * q[2] - q[1] * q[3]);
  rot[2][3] = 2.0 * (q[4] * q[3] + q[1] * q[2]);
  rot[3][3] = q[1]*q[1] - q[2]*q[2] - q[3]*q[3] + q[4]*q[4];
}

/* Rotate strucutre. This is to be done with all the atoms of the other strucutre.*/
void Rotate_strucutre(double *x2, double *y2, double *z2, int n, double rot[4][4])
{
  double xrot,yrot,zrot;
  int ia;

  // rotate the second structure
  for(ia=1; ia<=n; ia++)  
    {       
      xrot = x2[ia]*rot[1][1] + y2[ia]*rot[1][2] + z2[ia]*rot[1][3];
      yrot = x2[ia]*rot[2][1] + y2[ia]*rot[2][2] + z2[ia]*rot[2][3];
      zrot = x2[ia]*rot[3][1] + y2[ia]*rot[3][2] + z2[ia]*rot[3][3];
      x2[ia] = xrot;
      y2[ia] = yrot;
      z2[ia] = zrot;
    }
}

/* Translate back. This is to be done with all the atoms of both the strucutres*/
void Translate_back(double *xa, double *ya, double *za, double *xb, double *yb, double *zb, int AtomNum, double midxa, double midya, double midza)
{
  int i;
  
  for(i=1; i<=AtomNum; i++)       
    {
      xa[i]+=midxa;
      ya[i]+=midya;
      za[i]+=midza;
      
      xb[i]+=midxa;
      yb[i]+=midya;
      zb[i]+=midza;
    }
}


/* Align two strucutres using the coordinates of a subset of atoms. After alignment
 * the coordinates of the other strucutre is overwritten with new values.
 * 
 * Lei uses i-based arrays for the heart of the calcaultion. It is important to keep this
 * straight.
 *
 * Indices of atoms (of the subset that is to be used for alignment) are 1-based.
 * The array that stores these indices is a regular C-array with index starting from 0.*/
void Align_two_strucutres_using_subset_lei(FR_DAT *fr_ref, int Num_atoms_in_subset, int *Indices_atoms_in_subset, FR_DAT *fr)
{
  int i, Num_atoms;
  double *Reference_x, *Reference_y, *Reference_z;
  double *Other_x, *Other_y, *Other_z;
  double *Reference_x_subset, *Reference_y_subset, *Reference_z_subset;
  double *Other_x_subset, *Other_y_subset, *Other_z_subset;
  double midxa, midya, midza, midxb, midyb, midzb;
  double rot[4][4];

  Num_atoms = fr_ref->natoms;
  /* Memory allocation. All of these array have indices that start from 1. */
  Reference_x = (double*)malloc((Num_atoms+1) * sizeof(double));
  Reference_y = (double*)malloc((Num_atoms+1) * sizeof(double));
  Reference_z = (double*)malloc((Num_atoms+1) * sizeof(double));
  Other_x = (double*)malloc((Num_atoms+1) * sizeof(double));
  Other_y = (double*)malloc((Num_atoms+1) * sizeof(double));
  Other_z = (double*)malloc((Num_atoms+1) * sizeof(double));
  Reference_x_subset = (double*)malloc((Num_atoms_in_subset+1) * sizeof(double));
  Reference_y_subset = (double*)malloc((Num_atoms_in_subset+1) * sizeof(double));
  Reference_z_subset = (double*)malloc((Num_atoms_in_subset+1) * sizeof(double));
  Other_x_subset = (double*)malloc((Num_atoms_in_subset+1) * sizeof(double));
  Other_y_subset = (double*)malloc((Num_atoms_in_subset+1) * sizeof(double));
  Other_z_subset = (double*)malloc((Num_atoms_in_subset+1) * sizeof(double));

  /* Copy the coordiantes that are used for alighment into new data strucutres*/
  for(i = 0 ; i < Num_atoms ; i++)
    {
      Reference_x[i+1] = fr_ref->x[i][0];
      Reference_y[i+1] = fr_ref->x[i][1];
      Reference_z[i+1] = fr_ref->x[i][2];

      Other_x[i+1] = fr->x[i][0];
      Other_y[i+1] = fr->x[i][1];
      Other_z[i+1] = fr->x[i][2];
    }
  for(i = 0 ; i < Num_atoms_in_subset ; i++)
    {
      Reference_x_subset[i+1] = fr_ref->x[Indices_atoms_in_subset[i]-1][0];
      Reference_y_subset[i+1] = fr_ref->x[Indices_atoms_in_subset[i]-1][1];
      Reference_z_subset[i+1] = fr_ref->x[Indices_atoms_in_subset[i]-1][2];
      
      Other_x_subset[i+1] = fr->x[Indices_atoms_in_subset[i]-1][0];
      Other_y_subset[i+1] = fr->x[Indices_atoms_in_subset[i]-1][1];
      Other_z_subset[i+1] = fr->x[Indices_atoms_in_subset[i]-1][2];
    }
  
  /* Calculate the centers of mass of both the strucutres and move each of them.*/
  Calc_center(Reference_x_subset, Reference_y_subset, Reference_z_subset, Other_x_subset, Other_y_subset, Other_z_subset, Num_atoms_in_subset, &midxa, &midya, &midza, &midxb, &midyb, &midzb);
  /* Translate the subset*/
  Translate(Reference_x_subset, Reference_y_subset, Reference_z_subset, Other_x_subset, Other_y_subset, Other_z_subset, Num_atoms_in_subset, midxa, midya, midza, midxb, midyb, midzb);
  /* Translate all the atoms*/
  Translate(Reference_x, Reference_y, Reference_z, Other_x, Other_y, Other_z, Num_atoms, midxa, midya, midza, midxb, midyb, midzb);

  /* Calculate the rotation matrix*/
  Calc_rotation_matrix(Reference_x_subset, Reference_y_subset, Reference_z_subset, Other_x_subset, Other_y_subset, Other_z_subset, Num_atoms_in_subset, rot);
  
  /* Rotate and translate the other strucutre and translate the reference structure.*/
  /* Rotate*/
  Rotate_strucutre(Other_x, Other_y, Other_z, Num_atoms, rot);
  /* Translate back*/
  Translate_back(Reference_x, Reference_y, Reference_z, Other_x, Other_y, Other_z, Num_atoms, midxa, midya, midza);

  /* Copy new positions into input data strucutres*/
  for(i = 0 ; i < Num_atoms ; i++)
    {
      fr_ref->x[i][0] = Reference_x[i+1];
      fr_ref->x[i][1] = Reference_y[i+1];
      fr_ref->x[i][2] = Reference_z[i+1];

      fr->x[i][0] = Other_x[i+1];
      fr->x[i][1] = Other_y[i+1];
      fr->x[i][2] = Other_z[i+1];
    }


  free(Reference_x);
  free(Reference_y);
  free(Reference_z);
  free(Other_x);
  free(Other_y);
  free(Other_z);
  free(Reference_x_subset);
  free(Reference_y_subset);
  free(Reference_z_subset);
  free(Other_x_subset);
  free(Other_y_subset);
  free(Other_z_subset);
}




/* Read the indices in the subset of atoms used for alignemnt. The indices are 1-based.
 * Filename: 'SUBSET_FOR_ALIGNMENT'
 * Format:
 * ---------------------------------
 *  Number of atoms in the subset
 *  Atom indices one per line. Indices are 1-based meaning they start from 1.
 * ---------------------------------
 */
void Read_alloc_subset_for_alignment(int Num_atoms, int *Num_atoms_in_subset, int **Indices_atoms_in_subset)
{
  FILE *in;
  int lin_num, i;

  in = fopen("SUBSET_FOR_ALIGNMENT", "r");
  if(in == NULL)
    {
      printf("File 'SUBSET_FOR_ALIGNMENT' not found\n");
      exit(1);
    }
  lin_num = 1;
  if(fscanf(in, "%d", Num_atoms_in_subset) != 1)
    {
      printf("Error in line %d of file 'SUBSET_FOR_ALIGNMENT' \n", lin_num);
      exit(1);
    }
  if(*Num_atoms_in_subset > Num_atoms || *Num_atoms_in_subset < 1)
    {
      printf("Number of atoms in the subset can not be larger than the number of atoms in the system or it can not be smaller than 1.\n");
      exit(1);
    }
  lin_num++;
  /* Allocate memory*/
  *Indices_atoms_in_subset = (int*)malloc((*Num_atoms_in_subset) * sizeof(int));
  for(i = 0 ; i < (*Num_atoms_in_subset) ; i++)
    {
      if(fscanf(in, "%d", &(*Indices_atoms_in_subset)[i]) != 1)
        {
          printf("Error in line %d of file 'SUBSET_FOR_ALIGNMENT' \n", lin_num);
          exit(1);
        }
      if((*Indices_atoms_in_subset)[i] > Num_atoms || (*Indices_atoms_in_subset)[i] < 1)
        {
          printf("Invalid index in line %d of 'SUBSET_FOR_ALIGNMENT'\n", lin_num);
          exit(1);
        }
      lin_num++; 
    }
  fclose(in);  
}



