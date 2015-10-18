//  "Generalized-restart Lanczos Method for Schrödinger Equation in Semiconductor Simulation"
//  Created by Hantian Zhang, Xiaolin Guo on 10/09/15.
//  Copyright (c) 2015 Hantian Zhang, Xiaolin Guo. All rights reserved.

#ifndef UTIL
#define UTIL

#include <list>
#include <vector>

struct maRowEle {
	int col;
	double value;
	maRowEle(int r, double v): col(r), value(v) {};
};

struct CSRMatrix_t
{
	double *values;
	int    *colIdx;
	int    *rowStart;
};

struct CSRLocalMatrix_t
{
	CSRMatrix_t block;
	CSRMatrix_t ghost0;
	CSRMatrix_t ghost1;
};

struct SymmSeptDiagLocalMatrix_t
{
	double *block[4];
	double *ghost0[2];
	double *ghost1[2];
};

typedef std::vector< std::list<maRowEle> > triplet_m;

void mat_mul(double **A, double **B, double **C, int M, int N, int K);
double dis(double *A, double *B, int n);

void septdiagmvm(int M, int n, double* A[4], double* x, double* y);
void ghostmvm(int M, int n0, int n1, double* A0[2], double* A1[2],
              double* x0, double* x1, double* y);

void smvm(int m, const double* values, const int* col_idx, const int* row_start, double* x, double* y);
void smvm(int m, CSRMatrix_t& A, double* x, double* y);
void smvm_mpi(int                  N,
              int                       M,
              CSRLocalMatrix_t          &csr_loc,
              SymmSeptDiagLocalMatrix_t &ssd_loc,
              const double*             x,
              double*                   y);

void ToCSR(int N,  double* values,  int* col_idx, int* row_start, double **matrix);
void printCSR(int N,  double* values,  int* col_idx, int* row_start);

triplet_m CSRToMA(int N, double*& values, int*& colIdx, int*& rowStart);
void maToCSR(triplet_m& ma,
             double*& values, int*& colIdx, int*& rowStart);
void addToMA(double a, int row, int col, triplet_m& ma);

void lambdaAplusB(bool isPlus, int N, double lambda,
                  double *A_values, int *A_colIdx, int *A_rowStart,
                  double *B_values, int *B_colIdx, int *B_rowStart,
                  double*& C_values, int*& C_colIdx, int*& C_rowStart
                 );

double dotprod(const double* x, const double* y, int N);

double dotprodSIMD(const double* x, const double* y, int N);

void x_minusequal_beta_times_y(double* x, double beta, double* y, int N);

#endif
