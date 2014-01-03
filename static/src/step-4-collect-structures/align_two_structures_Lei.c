/* This code aligns tw ostructures in vmd xyz format.
 * It uses Lei Huang's code for alignment.*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/*This is from Lei's code.*/
#define ZERO	1.0E-10

#include "align_two_structures_Lei.h"

/******** These routines are from Lei Huang code 'superpose.cpp'. I made very few modifications. ***********/

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


void quatfit (double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n)
{
	double xxyx,xxyy,xxyz, xyyx,xyyy,xyyz, xzyx,xzyy,xzyz;
	double xrot,yrot,zrot, work1[5],work2[5], q[5],d[5], rot[4][4], c[5][5], v[5][5];
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

	for(ia=1; ia<=n; ia++)	{
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

	
	for(ia=1; ia<=n; ia++)	{	// rotate the second structure
		xrot = x2[ia]*rot[1][1] + y2[ia]*rot[1][2] + z2[ia]*rot[1][3];
		yrot = x2[ia]*rot[2][1] + y2[ia]*rot[2][2] + z2[ia]*rot[2][3];
		zrot = x2[ia]*rot[3][1] + y2[ia]*rot[3][2] + z2[ia]*rot[3][3];
		x2[ia] = xrot;
		y2[ia] = yrot;
		z2[ia] = zrot;
	}

	return;
}

/** COMMENTS BY AVISEK: I commented out the rmsd calculation so the routine does not return anything.
    (xa, ya, za) are arrays for the reference structure and (xb, yb, zb) are those for the other structure.
**/
// the index starting from 1 !!!
void CalRMSD(double *xa, double *ya, double *za, double *xb, double *yb, double *zb, int AtomNum)	
{
        //double RMSD;
	double midxa, midya, midza, midxb, midyb, midzb;
	int i;

	midxa=midya=midza=midxb=midyb=midzb=0.0;
	for(i=1; i<=AtomNum; i++)	{
		midxa += (xa[i]);
		midya += (ya[i]);
		midza += (za[i]);

		midxb += (xb[i]);
		midyb += (yb[i]);
		midzb += (zb[i]);
	}
	midxa/=AtomNum;
	midya/=AtomNum;
	midza/=AtomNum;
	midxb/=AtomNum;
	midyb/=AtomNum;
	midzb/=AtomNum;

	for(i=1; i<=AtomNum; i++)	{
		xa[i]-=midxa;
		ya[i]-=midya;
		za[i]-=midza;
		xb[i]-=midxb;
		yb[i]-=midyb;
		zb[i]-=midzb;
	}

	quatfit(xa, ya, za, xb, yb, zb, AtomNum);
	//RMSD=rmsfit(xa, ya, za, xb, yb, zb, AtomNum);


	for(i=1; i<=AtomNum; i++)	{	//restore structure1, and rotate + translate structure2 
		xa[i]+=midxa;
		ya[i]+=midya;
		za[i]+=midza;

		xb[i]+=midxa;
		yb[i]+=midya;
		zb[i]+=midza;
	}

	//return RMSD;
}

/**************************** End of Lei Huang's routines ******************************/




/* Align two data sets. The reference array(s) has dimension one more than the number of atoms.
 * This is necause Lei uses 1-based arrays, not 0-based.
 * THIS IS IMPORTANT TO KEEP STRAIGHT. LEI'S ROUTINE WORKS WITH ARRAYS WHOSE INDEX START FROM
 * 1 NOT ZERO.
 *
 * After alignemnt the coordinates are other structures are overwritten with new coordinates.
 * The reference arrays should have been allocated before.*/
void Align_two_structures_lei(double *Reference_x, double *Reference_y, double *Reference_z, FR_DAT *other_fr)
{
  int i;
  double *Other_x, *Other_y, *Other_z;

  /*Lei's routine uses arrays whose indices start from 1.*/
  Other_x = (double*)malloc((other_fr->natoms+1) * sizeof(double));
  Other_y = (double*)malloc((other_fr->natoms+1) * sizeof(double));
  Other_z = (double*)malloc((other_fr->natoms+1) * sizeof(double));

  /*Copy the coordinates of the other structure*/
  for(i = 0 ; i < other_fr->natoms ; i++)
    {
      Other_x[i+1] = other_fr->x[i][0];
      Other_y[i+1] = other_fr->x[i][1];
      Other_z[i+1] = other_fr->x[i][2]; 
    }

  /*This routine is from Lei and does the alignment*/
  CalRMSD(Reference_x, Reference_y, Reference_z, Other_x, Other_y, Other_z, other_fr->natoms);
  
  /*Copy new coordinates back*/
  for(i = 0 ; i < other_fr->natoms ; i++)
    {
      other_fr->x[i][0] = Other_x[i+1];
      other_fr->x[i][1] = Other_y[i+1];
      other_fr->x[i][2] = Other_z[i+1];
    }

  free(Other_x);
  free(Other_y);
  free(Other_z);
}









