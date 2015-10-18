//  "Generalized-restart Lanczos Method for Schrödinger Equation in Semiconductor Simulation"
//  "Section 4.3: Distributed-memory parallelism"
//  Created by Hantian Zhang, Xiaolin Guo on 10/09/15.
//  Copyright (c) 2015 Hantian Zhang, Xiaolin Guo. All rights reserved.

#include "CG_mpi.h"
#include "util.h"

#ifdef TIMER
#include "timer.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <memory.h>
#include <mpi.h>
#include "x86intrin.h"

using namespace std;
void get_csr_loc(int              N,
                 triplet_m        &ma_block,
                 triplet_m        &ma_ghost0,
                 triplet_m        &ma_ghost1,
                 CSRLocalMatrix_t &csr_loc);

void CG_mpi(CSRLocalMatrix_t          &csr_loc,
            double                    *&b_loc,
            double                    *&x_loc,
            int                       N,
            int                       M,
            int                       n,
            double                    threshold,
            int                       rank,
            int                       procs)
{
#ifdef TIMER
	double start = MPI_Wtime();
#endif

	MPI_Request req[4];
	MPI_Status status[4];

	int maxiter = N; // for now, need to be modified
	threshold = threshold * threshold; // we use the square of norm

	int d_loc_size(M + 2 * n); // size with ghost region included
	if (rank == 0 || rank == procs - 1) d_loc_size = M + n;
	if (procs == 1) d_loc_size = M;

	int d_offset(n); // where the non-ghost part starts
	if (rank == 0) d_offset = 0;

	for (int i = 0; i < d_loc_size; i++)
		x_loc[i] = 1. / N; // what is a better initial value?

	double *r_loc = (double*)malloc(sizeof(double) * M);
	double *d_loc = (double*)malloc(sizeof(double) * (d_loc_size));
	double *Ad = (double*)malloc(sizeof(double) * M); // A * d(:,k)

	double ghostvec0[d_offset], ghostvec1[d_loc_size - d_offset - M];
	smvm(M, csr_loc.block, x_loc + d_offset, r_loc);

	smvm(d_offset, csr_loc.ghost0, x_loc, ghostvec0);
	for (int i = 0; i < d_offset; i++)
		r_loc[i] += ghostvec0[i];

	smvm(d_loc_size - d_offset - M, csr_loc.ghost1, x_loc + d_offset + M, ghostvec1);
	for (int i = 0; i < d_loc_size - d_offset - M; i++)
		r_loc[M - n + i] += ghostvec1[i];

	for (int i = 0; i < M; i++)
		r_loc[i] = b_loc[i + d_offset] - r_loc[i]; // r(:,1) = b - A*x0

	if (rank > 0) {
		MPI_Irecv(d_loc, n, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &req[0] );
		MPI_Isend(r_loc, n, MPI_DOUBLE, rank - 1, 2, MPI_COMM_WORLD, &req[1] );
	}
	else {
		req[0] = MPI_REQUEST_NULL;
		req[1] = MPI_REQUEST_NULL;
	}

	if (rank < procs - 1) {
		MPI_Irecv(d_loc + d_offset + M, n, MPI_DOUBLE, rank + 1, 2, MPI_COMM_WORLD, &req[2] );
		MPI_Isend(r_loc + M - n, n, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &req[3] );
	}
	else {
		req[2] = MPI_REQUEST_NULL;
		req[3] = MPI_REQUEST_NULL;
	}

	memcpy(d_loc + d_offset, r_loc, sizeof(double) * M); // d(:,1) = r(:,1)

	// these two should only be updated on root, or not as now ??
	// 3 gather + 1 broadcast vs 2 Allgather
	double dotprodrr = 0.; // r(:,k)' * r(:,k)
	double dotprodrr1 = 0.;

	double dotprodrr1_loc = dotprod(r_loc, r_loc, M); // r(:,k+1)' * r(:,k+1)

	MPI_Allreduce(&dotprodrr1_loc, &dotprodrr1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	for (int k = 0; k < maxiter; k++)
	{
		CG_counter += 1;

		dotprodrr = dotprodrr1;

		smvm(M, csr_loc.block, d_loc + d_offset, Ad);

		MPI_Waitall(4, req, status);

		smvm(d_offset, csr_loc.ghost0, d_loc, ghostvec0);
		for (int i = 0; i < d_offset; i++)
			Ad[i] += ghostvec0[i];
		smvm(d_loc_size - d_offset - M, csr_loc.ghost1, d_loc + d_offset + M, ghostvec1);
		for (int i = 0; i < d_loc_size - d_offset - M; i++)
			Ad[M - n + i] += ghostvec1[i];

		double alpha = 0.;

		double alpha_loc = dotprod(d_loc + d_offset, Ad, M);
	
		MPI_Allreduce(&alpha_loc, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		alpha = dotprodrr / alpha;

		for (int i = 0; i < M; ++i)
		{
			x_loc[i + d_offset] += alpha * d_loc[i + d_offset];
			r_loc[i] -= alpha * Ad[i]; // r is now r(:,k+1)
		}

		dotprodrr1_loc = dotprod(r_loc, r_loc, M);

		MPI_Allreduce(&dotprodrr1_loc, &dotprodrr1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		if (dotprodrr1 < threshold) break;

		double beta = dotprodrr1 / dotprodrr;

		for (int i = 0; i < M; i++)
			d_loc[i + d_offset] = r_loc[i] + beta * d_loc[i + d_offset]; // d is now d(:,k+1)

		if (rank > 0) {
			MPI_Irecv(d_loc, n, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &req[0] );
			MPI_Isend(d_loc + d_offset, n, MPI_DOUBLE, rank - 1, 2, MPI_COMM_WORLD, &req[1] );
		}
		else {
			req[0] = MPI_REQUEST_NULL;
			req[1] = MPI_REQUEST_NULL;
		}

		if (rank < procs - 1) {
			MPI_Irecv(d_loc + d_offset + M, n, MPI_DOUBLE, rank + 1, 2, MPI_COMM_WORLD, &req[2] );
			MPI_Isend(d_loc + d_offset + M - n, n, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &req[3] );
		}
		else {
			req[2] = MPI_REQUEST_NULL;
			req[3] = MPI_REQUEST_NULL;
		}

	}

#ifdef TIMER
	timer_CG_mpi += MPI_Wtime() - start;
#endif

	free(r_loc);
	free(d_loc);
	free(Ad);
}

void CG_mpi_wrapper(int                       N,
                    CSRLocalMatrix_t          &csr_loc,
                    double                    *b,
                    double                    *x,
                    double                    threshold)
{
	int rank, procs;
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (N % procs != 0)
		cout << "N % procs != 0" << endl;
	int M = N / procs;
	int n = sqrt(N);
	if (M < n)
		cout << "M < n" << endl;

	int d_loc_size(M + 2 * n); // size with ghost region included
	if (rank == 0 || rank == procs - 1) d_loc_size = M + n;
	if (procs == 1) d_loc_size = M;
	int d_offset(n); // where the non-ghost part starts
	if (rank == 0) d_offset = 0;

	double *b_loc = (double*)malloc(sizeof(double) * d_loc_size);
	double *x_loc = (double*)malloc(sizeof(double) * d_loc_size);

	memcpy(b_loc, b + rank * M - d_offset, sizeof(double) * d_loc_size);

	CG_mpi(csr_loc,
	       b_loc, x_loc,
	       N, M, n,
	       threshold,
	       rank, procs);

	MPI_Allgather(x_loc + d_offset, M, MPI_DOUBLE,
	           x, M, MPI_DOUBLE, MPI_COMM_WORLD);

	free(b_loc);
	free(x_loc);
}

void get_local_matrix(double                    *values,
                      int                       *colIdx,
                      int                       *rowStart,
                      int                       N,
                      CSRLocalMatrix_t          &csr_loc)
{
	int rank, procs;
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (N % procs != 0)
		cout << "N % procs != 0" << endl;
	int M = N / procs;
	int n = sqrt(N);
	if (M < n)
		cout << "M < n" << endl;

	triplet_m ma = CSRToMA(N, values, colIdx, rowStart);

	triplet_m ma_block(M);
	for (int i = 0; i < M; i++)
	{
		int colstart = rank * M;
		int colend   = colstart + M;
		for (auto ele : ma[i + rank * M])
		{
			if (ele.col >= colstart && ele.col < colend)
			{
				ele.col -= colstart;
				ma_block[i].push_back(ele);
			}
		}
	}

	triplet_m ma_ghost0(n);
	if (rank > 0)
		for (int i = 0; i < n; i++)
		{
			int colstart = rank * M - n;
			int colend   = colstart + n;
			for (auto ele : ma[i + rank * M])
			{
				if (ele.col >= colstart && ele.col < colend) {
					ele.col -= colstart;
					ma_ghost0[i].push_back(ele);
				}
			}
		}

	triplet_m ma_ghost1(n);
	if (rank < procs - 1)
		for (int i = 0; i < n; i++)
		{
			int colstart = (rank + 1) * M;
			int colend   = colstart + n;
			for (auto ele : ma[i + (rank + 1) * M - n])
			{
				if (ele.col >= colstart && ele.col < colend) {
					ele.col -= colstart;
					ma_ghost1[i].push_back(ele);
				}
			}
		}

	get_csr_loc(N, ma_block, ma_ghost0, ma_ghost1, csr_loc);

}

void get_csr_loc(int              N,
                 triplet_m        &ma_block,
                 triplet_m        &ma_ghost0,
                 triplet_m        &ma_ghost1,
                 CSRLocalMatrix_t &csr_loc)
{
	int rank, procs;
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	maToCSR(ma_block, csr_loc.block.values, csr_loc.block.colIdx, csr_loc.block.rowStart);

	if (rank > 0)
	{
		maToCSR(ma_ghost0,
		        csr_loc.ghost0.values, csr_loc.ghost0.colIdx, csr_loc.ghost0.rowStart);
	}


	if (rank < procs - 1)
	{
		maToCSR(ma_ghost1,
		        csr_loc.ghost1.values, csr_loc.ghost1.colIdx, csr_loc.ghost1.rowStart);
	}

}

