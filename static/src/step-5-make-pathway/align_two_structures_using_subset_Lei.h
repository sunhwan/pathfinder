#ifndef _align_two_structures_using_subset_Lei_h
#define _align_two_structures_using_subset_Lei_h

#include "md_ld_common_d_dimension_namd.h"


void jacobi (int n, int np, double a[5][5], double d[5], double v[5][5], double b[5], double z[5]);

void Calc_center(double *xa, double *ya, double *za, double *xb, double *yb, double *zb, int AtomNum, double *midxa, double *midya, double *midza, double *midxb, double *midyb, double *midzb);

void Translate(double *xa, double *ya, double *za, double *xb, double *yb, double *zb, int AtomNum, double midxa, double midya, double midza, double midxb, double midyb, double midzb);

void Calc_rotation_matrix(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n, double rot[4][4]);

void Rotate_strucutre(double *x2, double *y2, double *z2, int n, double rot[4][4]);

void Translate_back(double *xa, double *ya, double *za, double *xb, double *yb, double *zb, int AtomNum, double midxa, double midya, double midza);

void Align_two_strucutres_using_subset_lei(FR_DAT *fr_ref, int Num_atoms_in_subset, int *Indices_atoms_in_subset, FR_DAT *fr);


void Read_alloc_subset_for_alignment(int Num_atoms, int *Num_atoms_in_subset, int **Indices_atoms_in_subset);

#endif
