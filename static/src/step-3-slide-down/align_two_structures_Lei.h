#ifndef _align_two_structures_Lei_h
#define _align_two_structures_Lei_h


#include "md_ld_common_d_dimension_namd.h"


void jacobi (int n, int np, double a[5][5], double d[5], double v[5][5], double b[5], double z[5]);

void quatfit (double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n);

void CalRMSD(double *xa, double *ya, double *za, double *xb, double *yb, double *zb, int AtomNum);

void Align_two_structures_lei(double *Reference_x, double *Reference_y, double *Reference_z, FR_DAT *other_fr);


#endif
