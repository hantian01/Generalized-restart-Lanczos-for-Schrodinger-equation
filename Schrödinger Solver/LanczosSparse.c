//  "Generalized-restart Lanczos Method for Schrödinger Equation in Semiconductor Simulation"
//  "Section 2.1: Generalized restart Lanczos algorithm for GHEP"
//  Created by Hantian Zhang, Xiaolin Guo on 10/09/15.
//  Copyright (c) 2015 Hantian Zhang, Xiaolin Guo. All rights reserved.

#include "LanczosSparse.h"
#include "stdlib.h"
#include "qr.h"
#include "util.h"
#include "CG_mpi.h"
#include "mpi.h"

#include <string.h>
#include <iostream>
#include <cmath>
#ifdef TIMER
    #include "timer.h"
#endif
using namespace std;

void LanczosSparse(double *A_values,
                   int    *A_colIdx,
                   int    *A_rowStart,
                   double *M_values,
                   int    *M_colIdx,
                   int    *M_rowStart,
                   int    N,
                   int    M_dim,
                   double **V,
                   double **W,
                   double **T_lanc,
                   int    qr_iter_num,
                   double *eigenvalues,
                   double **ritz,
                   double threshold
                  )
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

    double r[N], alpha[M_dim], beta[M_dim], w_init[N], q_init[N], q[N], u[N];

    for (i = 0; i < N; ++i) {
        w_init[i] = 0;
        q_init[i] = sqrt(1.0 / N);
        u[i] = 0;
    }

    // r = M q_init
    smvm(N, M_values, M_colIdx, M_rowStart, q_init, r);

    // beta_0 = (q_init,r)^{1/2}
    double beta_0 = sqrt(fabs(dotprod(q_init, r, N)));

    //++++

    k = 1;

    while (k <= M_dim) {
        counter += 1;
        for (j = k; j <= M_dim ; ++j) {
            if (j == 1) {
                for (h = 0; h < N; h++) {
                    W[j - 1][h] = r[h] / beta_0;
                    V[j - 1][h] = q_init[h] / beta_0;
                }
            }
            else {
                for (h = 0; h < N; h++) {
                    W[j - 1][h] = r[h] / beta[j - 2];
                    V[j - 1][h] = q[h] / beta[j - 2];
                }
            }

            //r = A v_j
            smvm(N, A_values, A_colIdx, A_rowStart, V[j - 1], r);

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
            // Using CG to solve M q = r for q;
            CG_mpi_wrapper(N, M_csr_loc,
                           r, q, threshold);

            beta[j - 1] = sqrt(fabs(dotprod(q, r, N)));
            // Assemble T
            T_lanc[j - 1][j - 1] = alpha[j - 1];
            T_lanc[j - 1][j] = beta[j - 1];
            T_lanc[j][j - 1] = beta[j - 1];

           if(beta[j-1] == 0){
               if (j > 1)
              qr(T_lanc,j-1,j-1,eigenvalues,ritz,qr_iter_num);
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
            eps = fabs(beta[M_dim - 1] * ritz[M_dim - 1][k - 1]);

            if (fabs(eps) < threshold) {
                k++;
// // ++++++++  Lock convergent eigenpairs
                //u_k = V_m * s_k;
                memset(u, 0, sizeof(double) * N);
                for (j = 0; j < M_dim; ++j)
                    x_minusequal_beta_times_y(u, -ritz[j][k - 1], V[j], N);

                //q_k = u_k/norm(u_k);
                sum = sqrt(dotprod(u, u, N));

                for (i = 0; i < N ; ++i) {
                    q[i] = u[i] / sum;
                }

                // r = M q
                smvm(N, M_values, M_colIdx, M_rowStart, q, r);

                // M-orthogonalization
                for (i = 0; i <= k - 2; i++) {
                    GS_counter += 1;
                    //y = (v_i , r)
                    y = dotprod(V[i], r, N);

                    //r = r - w_i y
                    x_minusequal_beta_times_y(r, y, W[i], N);
                }

                CG_mpi_wrapper(N, M_csr_loc,
                               r, q, threshold);
                // beta_{k-1} = (q , r)^{1/2}
                beta[k - 2] = sqrt(fabs(dotprod(q, r, N)));
//+++++ Lock eigenpair end;

                if (k > M_dim - 1) {
                    delete [] M_csr_loc.block.values;
                    delete [] M_csr_loc.block.colIdx;
                    delete [] M_csr_loc.block.rowStart;

                    if (rank > 0)
                    {
                        delete [] M_csr_loc.ghost0.values;
                        delete [] M_csr_loc.ghost0.colIdx;
                        delete [] M_csr_loc.ghost0.rowStart;
                    }
                    if (rank < procs - 1)
                    {
                        delete [] M_csr_loc.ghost1.values;
                        delete [] M_csr_loc.ghost1.colIdx;
                        delete [] M_csr_loc.ghost1.rowStart;
                    }
                    return;
                }

            }
//+++++ ELSE
            else {
                flag = 0;
                //u_k = V_m * s_k;
                memset(u, 0, sizeof(double) * N);
                for (j = 0; j < M_dim; ++j)
                    x_minusequal_beta_times_y(u, -ritz[j][k - 1], V[j], N);

                //q_k = u_k/norm(u_k);
                sum = sqrt(dotprod(u, u, N));

                if (k == 1) {
                    for (i = 0; i < N ; ++i) {
                        q_init[i] = u[i] / sum;
                    }
                }
                else {
                    for (i = 0; i < N ; ++i) {
                        q[i] = u[i] / sum;
                    }
                }

                // r = M q
                if (k == 1)
                    smvm(N, M_values, M_colIdx, M_rowStart, q_init, r);
                else
                    smvm(N, M_values, M_colIdx, M_rowStart, q, r);

                // M-orthogonalization
                for (i = 0; i <= k - 2; i++) {
                    GS_counter += 1;
                    //y = (v_i , r)
                    y = dotprod(V[i], r, N);
                    //r = r - w_i y
                    x_minusequal_beta_times_y(r, y, W[i], N);
                }

                CG_mpi_wrapper(N, M_csr_loc,
                               r, q, threshold);

                // beta_{k-1} = (q , r)^{1/2}
                if (k == 1)
                    beta_0 = sqrt(fabs(dotprod(q_init, r, N)));
                else
                    beta[k - 2] = sqrt(fabs(dotprod(q, r, N)));


            }
        } // end while(flag == 1)

    }

}
