//  "Generalized-restart Lanczos Method for Schrödinger Equation in Semiconductor Simulation"
//  Created by Hantian Zhang, Xiaolin Guo on 10/09/15.
//  Copyright (c) 2015 Hantian Zhang, Xiaolin Guo. All rights reserved.

#include "math.h"
#include "stdlib.h"
#include "stdio.h"

void mod_grams(double **T, double **Q, double **R,int M,int N);
void qr(double ** A, int M, int N, double * eigenvalues, double **eigenvectors,int niters){
    int iter;
    int i,j,k;
    double **Q;
    double **R;
    double **U;
    double **T;
    double **B;
    double **P;
    
    Q = (double **)malloc(sizeof(double *) * M);
    R = (double **)malloc(sizeof(double *) * M);
    T = (double **)malloc(sizeof(double *) * M);
    U = (double **)malloc(sizeof(double *) * M);
    B = (double **)malloc(sizeof(double *) * M);
    P = (double **)malloc(sizeof(double *) * M);
    for(i = 0;i < M; ++i){
        
        Q[i] = (double *)malloc(sizeof(double) * M);
        R[i] = (double *)malloc(sizeof(double) * M);
        T[i] = (double *)malloc(sizeof(double) * M);
        U[i] = (double *)malloc(sizeof(double) * M);
        B[i] = (double *)malloc(sizeof(double) * M);

        P[i] = (double *)malloc(sizeof(double) * M);
        
        for(j = 0; j < N; ++j){
            Q[i][j] = 0;
            R[i][j] = 0;
            U[i][j] = 0;
            if(i == j){
                U[i][j] = 1;
            }
            //printf("%lf ",A[i][j]);
            T[i][j] = A[i][j];
        }
    }
    
    
    for(iter = 0; iter < niters; ++iter){
        //[Q,R]=mod_grams(T);
        
        for(i = 0;i < M;++i){
            for(j = 0;j <N; ++j){
                B[i][j] = T[i][j];
            }
        }
    
        
        mod_grams(B,Q,R,M,N);
        
        //T=R*Q;
        double sum = 0;
        for(i = 0; i < M; i++){
            for(j = 0; j < N; j++){
                for(k = 0; k<N; k++){
                    sum += R[i][k] * Q[k][j];
                }
                T[i][j] = sum;
                sum = 0;
            }
        }
        //U = U*Q

        sum = 0;
        for(i = 0; i < M; i++){
            for(j = 0; j < N; j++){
                for(k = 0; k<N; k++){
                    sum += U[i][k] * Q[k][j];
                }
                P[i][j] = sum;
                sum = 0;
            }
        }
        sum = 0;
        for(i = 0; i < M; i++){
            for(j = 0; j < N; j++){
                U[i][j] = P[i][j];
            }
        }
           
    }
    for(i = 0; i < M; ++i){
        eigenvalues[i] = T[i][i];
    }

    for(i = 0; i < M; ++i){
        for(j = 0; j < N ; ++j){
            eigenvectors[i][j] = U[i][j];
        }
    }

    for (i = 0; i < M; ++i)
    {
        free(Q[i]);
        free(R[i]);
        free(T[i]);
        free(U[i]);
        free(B[i]);
        free(P[i]);
    }
    free(Q);
    free(R);
    free(T);
    free(U);
    free(B);
    free(P);
}
void mod_grams(double **T, double **Q, double **R,int M,int N){
    int i,j,k,m;
    double tmp_norm;
    for(i = 0; i < M; i++){
        for(j = 0; j < N; j++){
            Q[i][j] = 0;
            R[i][j] = 0;
        }
        Q[i][0] = T[i][0];
    }
    R[0][0] = 1;
    for(k = 0; k < N; k++){
        //R(k,k) = norm(T(1:m,k));
        tmp_norm = 0;
        for(m = 0;m < M; m++){
            tmp_norm += T[m][k] * T[m][k];
        }
        R[k][k] = sqrt(tmp_norm);
        
        //Q(1:m,k) = T(1:m,k)/R(k,k);
        for(m = 0; m < M; m++){
            Q[m][k] = T[m][k] / R[k][k];
        }
        int mm = 0;
        for( j = k + 1; j < N; j++){
            //R(k,j) = Q(1:m,k)'*T(1:m,j);
            double temp_r = 0;
            for(m = 0; m < M; ++m)
            {
                temp_r += Q[m][k] * T[m][j];
            }
            R[k][j] = temp_r;
            //T(1:m,j) = T(1:m,j)-Q(1:m,k)*R(k,j);
            for(mm = 0; mm < M; ++mm){
                T[mm][j] = T[mm][j] - Q[mm][k]*R[k][j];
            }
        }
    }
    
}
