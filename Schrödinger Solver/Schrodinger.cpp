//  "Generalized-restart Lanczos Method for Schrödinger Equation in Semiconductor Simulation"
//  "Section 3: Extraction of m lowest eigenstates"
//  Created by Hantian Zhang, Xiaolin Guo on 10/09/15.
//  Copyright (c) 2015 Hantian Zhang, Xiaolin Guo. All rights reserved.

#include "LanczosSparse.h"
#include "LanczosInverseSparse.h"
#include "Sch_Matrix.h"
#include "util.h"
#include "timer.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

#include <mpi.h>

using std::cout;
using std::endl;


void Schrodinger(int N_Sch, int Nx_Sch, int Ny_Sch,
                      int M_dim, int qr_iter,
                      double THRESHOLD, double THRESHOLD_inverse,
                      double *Domain)
{
    int id = MPI::COMM_WORLD.Get_rank();

    double sum;
    int i, j, k;

    counter = 0;
    InvCounter = 0;
    GS_counter = 0;
    CG_counter = 0;

    // Domain declaration
    double x_a, x_b, y_a, y_b;
    x_a = Domain[0]; x_b = Domain[1];
    y_a = Domain[2]; y_b = Domain[3];

    int Nx_init = 64;
    int Ny_init = 64;
    int N_Init = Nx_init*Ny_init;

    double eigen_value[M_dim];
    double *eigen_vector[M_dim];
    for (int i = 0; i < M_dim; ++i)
        eigen_vector[i] = new double [N_Sch];

    // Initialization
    double *A_values,  *M_values, *A_tildeVal;
    int *A_colIdx, *A_tildecolIdx, *M_colIdx;
    int *A_rowStart, *A_tilderowStart, *M_rowStart;

    double **ritz_vector = (double **)malloc(sizeof(double *) * M_dim);
    double **T_lanc = (double **)malloc(sizeof(double*) * (M_dim + 1) );
    double **V = (double **)malloc(sizeof(double*) * (M_dim));
    double **W = (double **)malloc(sizeof(double*) * (M_dim));
    for (i = 0; i < M_dim; ++i)
    {
        V[i] = (double *)malloc(sizeof(double) * N_Sch);
        W[i] = (double *)malloc(sizeof(double) * N_Sch);
        ritz_vector[i] = (double *)malloc(sizeof(double) * M_dim);
    }

    double *error = (double *)malloc(sizeof(double) * M_dim);
    double *Ax = (double *)malloc(sizeof(double) * N_Sch);
    double *Mx = (double *)malloc(sizeof(double) * N_Sch);
    double *diff = (double *)malloc(sizeof(double) * N_Sch);
    double *eigen_vecArr = (double *)malloc(sizeof(double) * N_Sch);

    for (i = 0; i < (M_dim + 1); i++) {
        T_lanc[i] = (double *)malloc(sizeof(double) * (M_dim + 1) );
    }

    GalerkinMatrix_Sch(A_values, A_colIdx, A_rowStart,
                       M_values, M_colIdx, M_rowStart,
                       Nx_init, Ny_init, x_a, x_b, y_a, y_b); // Compute A matrix

    counter++;
    if (id == 0)
    {
        // printf("Working counter = %d \n", counter);
        printf("Mesh nodes of Schrodinger = %d \n", N_Sch);
    }

// // +++++++++++++++++++++++++++++++++=
    // Using Lanczos to extract the extreme eigenvalue (the one has largest absolute value)
    LanczosSparse(A_values, A_colIdx, A_rowStart,
                  M_values, M_colIdx, M_rowStart,
                  N_Init, 2, V, W, T_lanc, qr_iter,
                  eigen_value, ritz_vector, THRESHOLD);
    // compute upper bound
    double tau = fabs(eigen_value[0]) + 1;
    // printf("Upper Bound (Tau) = %f \n", tau);
    // compute Eigen-vector corresponds to tau
    for (i = 0; i < N_Init; ++i) {
        sum = 0;
        for (j = 0; j < 2; ++j) {
            sum += V[j][i] * ritz_vector[j][0];
        }
        eigen_vector[0][i] = sum;
    }
    // Compute Error Norm
    for (i = 0; i < N_Init; ++i) {
        eigen_vecArr[i] = eigen_vector[0][i];
    }
    smvm(N_Init, A_values, A_colIdx, A_rowStart, eigen_vecArr, Ax);
    smvm(N_Init, M_values, M_colIdx, M_rowStart, eigen_vecArr, Mx);
    for (i = 0; i < N_Init; ++i) {
        diff[i] = Ax[i] - eigen_value[0] * Mx[i];
    }
    // sum = 0;
    // for (i = 0; i < N_Init; ++i) {
    //     sum += diff[i] * diff[i];
    // }
    // double error_tau = sqrt(sum);
    // printf("Error of Tau = %f \n", error_tau);
    // getchar();
//++++++++++++++++++++++++++++++++++++++
    // A_tilde = tau M - A , ensure the positive definite property of A_tilde
    lambdaAplusB(false, N_Init, tau,
                 M_values, M_colIdx, M_rowStart,
                 A_values, A_colIdx, A_rowStart,
                 A_tildeVal, A_tildecolIdx, A_tilderowStart
                );

    //Compute lowest eigenvalue
    LanczosSparse(A_tildeVal, A_tildecolIdx, A_tilderowStart,
                  M_values, M_colIdx, M_rowStart,
                  N_Init, 2, V, W, T_lanc, qr_iter,
                  eigen_value, ritz_vector, THRESHOLD);

    double LowEigenVal = tau - eigen_value[0];
    if (id == 0){
    printf("Lowest eigenvalue = %f \n", LowEigenVal );
    }  
    // compute Eigen-vector corresponds to tau
    for (i = 0; i < N_Init; ++i) {
        sum = 0;
        for (j = 0; j < 2; ++j) {
            sum += V[j][i] * ritz_vector[j][0];
        }
        eigen_vector[0][i] = sum;
    }
    // Compute Error Norm
    for (i = 0; i < N_Init; ++i) {
        eigen_vecArr[i] = eigen_vector[0][i];
    }
    smvm(N_Init, A_tildeVal, A_colIdx, A_rowStart, eigen_vecArr, Ax);
    smvm(N_Init, M_values, M_colIdx, M_rowStart, eigen_vecArr, Mx);
    for (i = 0; i < N_Init; ++i) {
        diff[i] = Ax[i] - eigen_value[0] * Mx[i];
    }
    // sum = 0;
    // for (i = 0; i < N_Init; ++i) {
    //     sum += diff[i] * diff[i];
    // }
    // double error_LowEigen = sqrt(sum);
    // printf("Error of Lowest eigenvalue = %f \n", error_LowEigen);
    // getchar();
//++++++++++++++++++++++++++++++++++++++++
    GalerkinMatrix_Sch(A_values, A_colIdx, A_rowStart,
                       M_values, M_colIdx, M_rowStart,
                       Nx_Sch, Ny_Sch, x_a, x_b, y_a, y_b); // Compute A matrix
    // Minimum shift to ensure positive definie property
    double shift = 0.0;
    if (LowEigenVal < 0.0){
        shift = 2 * fabs(LowEigenVal);
    }
    // A_tilde = shift M + A , eigenvalue in form of (shift + lambda)
    delete [] A_tildeVal;
    delete [] A_tildecolIdx;
    delete [] A_tilderowStart;
    lambdaAplusB(true, N_Sch, shift,
                 M_values, M_colIdx, M_rowStart,
                 A_values, A_colIdx, A_rowStart,
                 A_tildeVal, A_tildecolIdx, A_tilderowStart
                );

// ++++++++++++++++++++++++++++++++++++++++++
    // Solving Schrodinger equation, inverse iteration to extract low eigenpairs
    // explicity restarted Inverse lanczos for GHEP, (M_dim - 1) pairs of eigenvalue/vector are valid
    // do not use the last pair for computing electric density
    LanczosInverseSparse(A_tildeVal, A_tildecolIdx, A_tilderowStart,
                         M_values, M_colIdx, M_rowStart,
                         N_Sch, M_dim, V, W, T_lanc, qr_iter, eigen_value, ritz_vector, THRESHOLD_inverse);
    // compute Eigen-vector
    for (k = 0; k < M_dim; k++) {
        for (i = 0; i < N_Sch; ++i) {
            sum = 0;
            for (j = 0; j < M_dim; ++j) {
                sum += W[j][i] * ritz_vector[j][k];
            }
            eigen_vector[k][i] = sum;
        }
    }
    // Compute Error Norm
    for (k = 0; k < M_dim; k++) {
        for (i = 0; i < N_Sch; ++i) {
            eigen_vecArr[i] = eigen_vector[k][i];
        }
        smvm(N_Sch, A_tildeVal, A_colIdx, A_rowStart, eigen_vecArr, Ax);
        smvm(N_Sch, M_values, M_colIdx, M_rowStart, eigen_vecArr, Mx);
        for (i = 0; i < N_Sch; ++i) {
            diff[i] = Ax[i] - 1.0 / eigen_value[k] * Mx[i];
        }
        sum = 0;
        for (i = 0; i < N_Sch; ++i) {
            sum += diff[i] * diff[i];
        }
        error[k] = sqrt(sum);
    }
    // convert eigen_value to lambda = tau - lambda_B
    for (i = 0; i < M_dim; ++i) {
        eigen_value[i] = -shift + 1.0 / eigen_value[i];
    }

    // print eigenvalue and ritz vector
    if (id == 0) {
        cout << endl << "eigenvalues: " << endl;
        for (i = 0; i < M_dim - 1; ++i) {
            printf("  Eigen_value %lf\n", eigen_value[i]);
            printf("  Error = %lf\n", error[i]);
        }
        cout << endl;
    }

    if (id == 0)
    {
        printf("\n");
        printf("Lanczos Counter = %d\n", counter );
        printf("Inv-Lanczos Counter = %d\n", InvCounter );
        printf("GS othorgonal counter = %d\n",GS_counter );
        printf("CG LSE Counter = %d \n", CG_counter);
        printf("\n");
    }
}
