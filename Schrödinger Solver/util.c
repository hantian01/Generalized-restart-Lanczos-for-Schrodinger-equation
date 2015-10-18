//  "Generalized-restart Lanczos Method for Schrödinger Equation in Semiconductor Simulation"
//  Created by Hantian Zhang, Xiaolin Guo on 10/09/15.
//  Copyright (c) 2015 Hantian Zhang, Xiaolin Guo. All rights reserved.

#include "util.h"

#ifdef TIMER
#include "timer.h"
#endif

#include <algorithm>
#include <math.h>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <mpi.h>
#include <string.h>

using namespace std;

void mat_mul(double **A, double **B, double **C, int M, int K, int N) {
	//M * K   K * N
	int m, n, k;
	for (m = 0; m < M; m++) {
		for (n = 0; n < N; n++) {
			C[m][n] = 0;
			for (k = 0; k < K; k++) {
				C[m][n] += A[m][k] * B[k][n];
			}
		}
	}
	return;
}

double dis(double *A, double *B, int n)
{
	int i;
	int res = 0;
	for (i = 0; i < n; i++)
		res += (A[i] - B[i]) * (A[i] - B[i]);
	return sqrt(res);
}

void smvm(int m, const double* values, const int* col_idx,
          const int* row_start, double* x, double* y)
{
#ifdef TIMER
	double start = MPI_Wtime();
#endif

	for (int i = 0; i < m; i++) {

		double d = 0; /* scalar replacement since reused */
		/* loop over non-zero elements in row i */
		for (int j = row_start[i]; j < row_start[i + 1]; j++){
			d += values[j] * x[col_idx[j]];
		}
		y[i] = d;
	}

#ifdef TIMER
	timer_smvm += MPI_Wtime() - start;
#endif

}

void smvm(int m, CSRMatrix_t& A, double* x, double* y)
{
	smvm(m, A.values, A.colIdx, A.rowStart, x, y);
}


void smvm_mpi(int                  N,
              int                       M,
              CSRLocalMatrix_t          &csr_loc,
              SymmSeptDiagLocalMatrix_t &ssd_loc,
              const double*             x,
              double*                   y)
{
	int rank, procs;
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int n = sqrt(N);
	int d_loc_size(M + 2 * n); // size with ghost region included
	if (rank == 0 || rank == procs - 1) d_loc_size = M + n;
	if (procs == 1) d_loc_size = M;
	int d_offset(n); // where the non-ghost part starts
	if (rank == 0) d_offset = 0;

	double *y_loc = (double*)_mm_malloc(sizeof(double) * M, 64 );
	double *x_loc = (double*)_mm_malloc(sizeof(double) * d_loc_size, 64 );
	memcpy(x_loc, x + rank * M - d_offset, sizeof(double) * d_loc_size);

	double ghostvec0[d_offset], ghostvec1[d_loc_size - d_offset - M];
	smvm(M, csr_loc.block, x_loc + d_offset, y_loc);
	smvm(d_offset, csr_loc.ghost0, x_loc, ghostvec0);
	for (int i = 0; i < d_offset; i++)
		y_loc[i] += ghostvec0[i];
	smvm(d_loc_size - d_offset - M, csr_loc.ghost1, x_loc + d_offset + M, ghostvec1);
	for (int i = 0; i < d_loc_size - d_offset - M; i++)
		y_loc[M - n + i] += ghostvec1[i];

	MPI_Allgather(y_loc, M, MPI_DOUBLE,
	              y, M, MPI_DOUBLE, MPI_COMM_WORLD);
}

void ToCSR(int N,  double * values,  int* col_idx, int* row_start, double **matrix) {
	int i, j;
	int counter = 0;

	for (i = 0; i < N; ++i)
	{
		row_start[i] = counter;
		for (j = 0; j < N; ++j)
		{
			if (matrix[i][j] != 0) {
				values[counter] = matrix[i][j];
				col_idx[counter] = j;
				counter += 1;
			}
		}
	}
	row_start[N] = counter;
}

// is it possible that there are some nonzeros in ma?
void maToCSR(triplet_m & ma,
             double*& values, int*& colIdx, int*& rowStart)
{

	//sort entries in each column by their row
	triplet_m::iterator it;
	for (auto& l : ma)
		l.sort([](const maRowEle & a, const maRowEle & b)
	{ return a.col < b.col; }
	      );

	rowStart = new int [ma.size() + 1];
	rowStart[0] = 0;
	int nnz(0), k(0);
	for ( it = ma.begin(), k = 1;
	        it != ma.end(); it++, k++) {
		nnz += it -> size();
		rowStart[k] = nnz;
	}

	colIdx = new int [nnz];
	values = new double [nnz];
	for (it = ma.begin(), k = 0; it != ma.end(); it++)
		for (std::list<maRowEle>::iterator it1 = it -> begin();
		        it1 != it -> end(); it1++, k++)
		{
			colIdx[k] = it1 -> col;
			values[k] = it1 -> value;
		}

}

void printCSR(int N, double * values, int* col_idx, int* row_start)
{
	using namespace std;
	cout << "values" << endl;
	for (int i = 0; i < row_start[N] - 1; i++)
		cout << values[i] << ", ";
	cout << values[row_start[N] - 1] << endl;

	cout << "colIdx" << endl;
	for (int i = 0; i < row_start[N] - 1; i++)
		cout << col_idx[i] << ", ";
	cout << col_idx[row_start[N] - 1] << endl;

	cout << "rowStart" << endl;
	for (int i = 0; i < N; i++)
		cout << row_start[i] << ", ";
	cout << row_start[N] << endl;
}

triplet_m CSRToMA(int N, double*& values, int*& colIdx, int*& rowStart)
{
	int nnz = rowStart[N];
	triplet_m ma(N);
	for (int k = 0, row = 0; k < nnz; k++)
	{
		while (k >= rowStart[row] && row <= N)
			row++;

		double a = values[k];
		double col = colIdx[k];

		ma[row - 1].push_back(maRowEle(col, a));
	}

	return ma;
}

void addToMA(double a, int row, int col, triplet_m & ma)
{
	if (a == 0) return;

	list<maRowEle>::iterator it =
	    find_if(ma[row].begin(), ma[row].end(),
	[&col] (const maRowEle & ele) { return ele.col == col; } );
	// [&item](const T&v) { return item == v; }
	//     find(ma[row].begin(), ma[row].end(), col);

	if (it != ma[row].end())
	{
		it -> value += a;
	}
	else
	{
		ma[row].push_back(maRowEle(col, a));
	}

}

void lambdaAplusB(bool isPlus, int N, double lambda,
                  double * A_values, int *A_colIdx, int *A_rowStart,
                  double * B_values, int *B_colIdx, int *B_rowStart,
                  double*& C_values, int*& C_colIdx, int*& C_rowStart
                 )
{
	double *adj_B_values = B_values;
	if (!isPlus)
	{
		int nnz = B_rowStart[N];
		adj_B_values = new double [nnz];
		for (int i = 0; i < nnz; i++)
			adj_B_values[i] = -B_values[i];
	}

	triplet_m ma =
	    CSRToMA(N, A_values, A_colIdx, A_rowStart);
	triplet_m mb =
	    CSRToMA(N, adj_B_values, B_colIdx, B_rowStart);
	triplet_m mc(N);

	for (int row = 0; row < N; row++)
	{
		list<maRowEle>& lc = mc[row];

		for (const maRowEle& mre_a : ma[row])
			lc.push_back(
			    maRowEle(mre_a.col, lambda * mre_a.value)
			);

		for (const maRowEle& mre_b : mb[row])
			addToMA(mre_b.value, row, mre_b.col, mc);

	}

	maToCSR(mc, C_values, C_colIdx, C_rowStart);

	if (!isPlus) delete [] adj_B_values;
}

double dotprod(const double* x, const double* y, int N)
{
	double sum = 0.;
	for (int i = 0; i < N; i++) {
		sum += x[i] * y[i];
	}
	return sum;

}

void x_minusequal_beta_times_y(double* x, double beta, double* y, int N)
{
	for (int h = 0; h < N; h++) {
		x[h] -= y[h] * beta;
	}

}

