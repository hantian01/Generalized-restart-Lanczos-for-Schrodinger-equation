//  "Generalized-restart Lanczos Method for Schrödinger Equation in Semiconductor Simulation"
//  "Inverse variant of generalized restart Lanczos algorithm"
//  Created by Hantian Zhang, Xiaolin Guo on 10/09/15.
//  Copyright (c) 2015 Hantian Zhang, Xiaolin Guo. All rights reserved.

#include "LanczosInverseSparse.h"
#include "stdlib.h"
#include "math.h"
#include "qr.h"
#include "util.h"
#include <mpi.h>
#include <string.h>
#include <iostream>
#include "CG_mpi.h"
#ifdef TIMER
    #include "timer.h"
#endif

using namespace std;

void LanczosInverseSparse(double *A_values, int *A_colIdx, int *A_rowStart,
                          double *M_values, int *M_colIdx, int *M_rowStart,
                          int N, int M_dim, double **V, double **W,
                          double **T_lanc, int qr_iter_num, double *eigenvalues,
                          double **ritz, double threshold)
{

    CSRLocalMatrix_t A_csr_loc, M_csr_loc;
    get_local_matrix(A_values, A_colIdx, A_rowStart, N,
                     A_csr_loc);
    get_local_matrix(M_values, M_colIdx, M_rowStart, N,
                     M_csr_loc);
    int rank, procs;
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int i, j, k;
    int flag;
    int h;
    double sum, y;

    for (i = 0; i < M_dim; ++i) {
        memset(ritz[i], 0, sizeof(double) * M_dim);
        memset(T_lanc[i], 0, sizeof(double) * M_dim);
    }

    // double r[N], alpha[M_dim], beta[M_dim], w_init[N], r_init[N], q[N],  u[N];
    double *r = (double *)malloc(sizeof(double) * N);
    double *alpha = (double *)malloc(sizeof(double) * M_dim);
    double *beta = (double *)malloc(sizeof(double) * M_dim);
    double *w_init = (double *)malloc(sizeof(double) * N);
    double *r_init = (double *)malloc(sizeof(double) * N);
    double *q = (double *)malloc(sizeof(double) * N);

    double *u = (double *)malloc(sizeof(double ) * N);
    for (i = 0; i < N; ++i) {
        w_init[i] = 0;
        r_init[i] = sqrt(1.0 / N);
        u[i] = 0;
    }
    // q = M r_init

    smvm(N, M_values, M_colIdx, M_rowStart, r_init, q);

    // beta_0 = (q_init,r)^{1/2}
    double beta_0 = sqrt(fabs(dotprod(r_init, q, N)));


    k = 1;

    while (k <= M_dim) {
        InvCounter += 1;
        for (j = k; j <= M_dim ; ++j) {
            if (j == 1) {
                for (h = 0; h < N; h++) {
                    W[j - 1][h] = r_init[h] / beta_0;
                    V[j - 1][h] = q[h] / beta_0;
                }
            }
            else {
                for (h = 0; h < N; h++) {
                    W[j - 1][h] = r[h] / beta[j - 2];
                    V[j - 1][h] = q[h] / beta[j - 2];
                }
            }

            // Using CG to solve A r = v_j for r
            // printf("before LanzLinear for inverse iteration\n");
            CG_mpi_wrapper(N, A_csr_loc,
                           V[j - 1], r, threshold);

            // r =r - w_{j-1} beta_{j-1}
            if (j == 1)
                x_minusequal_beta_times_y(r, beta_0, w_init, N);
            else
                x_minusequal_beta_times_y(r, beta[j - 2], W[j - 2], N);

            // alpha_j = (v_j , r);
            alpha[j - 1] = dotprod(V[j - 1], r, N);
            // r = r - w_j alpha_j
            if (j == 1)
                x_minusequal_beta_times_y(r, alpha[j - 1], w_init, N);
            else
                x_minusequal_beta_times_y(r, alpha[j - 1], W[j - 1], N);

            // Gram-Schmidt orthogonalization
            for (i = 1; i <= j; ++i) {
                GS_counter += 1;
                //y = (v_i , r);
                y = dotprod(V[i - 1], r, N);

                //r = r - w_i y
                x_minusequal_beta_times_y(r, y, W[i - 1], N);
            }

            // q = M r;
            smvm(N, M_values, M_colIdx, M_rowStart, r, q);

            beta[j - 1] = sqrt(fabs(dotprod(q, r, N)));
            // Assemble T
            T_lanc[j - 1][j - 1] = alpha[j - 1];
            T_lanc[j - 1][j] = beta[j - 1];
            T_lanc[j][j - 1] = beta[j - 1];

            if (beta[j - 1] == 0) {
                qr(T_lanc, j - 1, j - 1, eigenvalues, ritz, qr_iter_num);
                return ;
            }//endif
        }//endfor


        //solve m eigenvalues.
        qr( T_lanc , M_dim , M_dim , eigenvalues, ritz, qr_iter_num);
        flag = 1;

        while (flag == 1) {
//++++++ Ckech residual norm
            double eps = 0;
            // eps = hm+1m s(m)k;
            eps = fabs(beta[M_dim - 1] * ritz[M_dim - 1][k - 1] / eigenvalues[k - 1]);

            if (fabs(eps) < threshold) {
                // printf("YES\n");
                k++;
// // ++++++++  Lock convergent eigenpairs
                //u_k = W_m * s_k;
                memset(u, 0, sizeof(double) * N);
                for (j = 0; j < M_dim; ++j)
                    x_minusequal_beta_times_y(u, -ritz[j][k - 1], W[j], N);

                //r_k = u_k/norm(u_k);
                sum = sqrt(dotprod(u, u, N));

                for (i = 0; i < N ; ++i) {
                    r[i] = u[i] / sum;
                }

                // M-orthogonalization
                for (i = 0; i <= k - 2; i++) {
                    GS_counter += 1;
                    //y = (v_i , r)
                    y = dotprod(V[i], r, N);

                    //r = r - w_i y
                    x_minusequal_beta_times_y(r, y, W[i], N);
                }

                // q = M r
                smvm(N, M_values, M_colIdx, M_rowStart, r, q);

                // beta_{k-1} = (q , r)^{1/2}
                beta[k - 2] = sqrt(fabs(dotprod(q, r, N)));
//+++++ Lock eigenpair end;

                if (k > M_dim - 1) {
                            free(r);
                            free(alpha);
                            free(beta);
                            free(w_init);
                            free(r_init);
                            free(q);
                            free(u);

                            delete [] A_csr_loc.block.values;
                            delete [] A_csr_loc.block.colIdx;
                            delete [] A_csr_loc.block.rowStart;
                            if (rank > 0)
                            {
                                delete [] A_csr_loc.ghost0.values;
                                delete [] A_csr_loc.ghost0.colIdx;
                                delete [] A_csr_loc.ghost0.rowStart;
                            }
                            if (rank < procs - 1)
                            {
                                delete [] A_csr_loc.ghost1.values;
                                delete [] A_csr_loc.ghost1.colIdx;
                                delete [] A_csr_loc.ghost1.rowStart;
                            }
                    return;
                }

            }
//+++++ ELSE
            else {
                flag = 0;
                memset(u, 0, sizeof(double) * N);
                //u_k = W_m * s_k;
                for (j = 0; j < M_dim; ++j)
                    x_minusequal_beta_times_y(u, -ritz[j][k - 1], W[j], N);

                //q_k = u_k/norm(u_k);
                sum = sqrt(dotprod(u, u, N));

                if (k == 1) {
                    for (i = 0; i < N ; ++i) {
                        r_init[i] = u[i] / sum;
                    }
                }
                else {
                    for (i = 0; i < N ; ++i) {
                        r[i] = u[i] / sum;
                    }
                }


                // M-orthogonalization
                for (i = 0; i <= k - 2; i++) {
                    GS_counter += 1;
                    //y = (v_i , r)
                    y = dotprod(V[i], r, N);
                    //r = r - w_i y
                    x_minusequal_beta_times_y(r, y, W[i], N);
                }

                // q = M r
                if (k == 1)
                    smvm(N, M_values, M_colIdx, M_rowStart, r_init, q);
                else
                    smvm(N, M_values, M_colIdx, M_rowStart, r, q);

                // beta_{k-1} = (q , r)^{1/2}
                if (k == 1)
                    beta_0 = sqrt(fabs(dotprod(r_init, q, N)));
                else
                    beta[k - 2] = sqrt(fabs(dotprod(r, q, N)));

            }
        }
    }

}
