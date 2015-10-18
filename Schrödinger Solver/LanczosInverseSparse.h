//  "Generalized-restart Lanczos Method for Schrödinger Equation in Semiconductor Simulation"
//  "Inverse variant of generalized restart Lanczos algorithm"
//  Created by Hantian Zhang, Xiaolin Guo on 10/09/15.
//  Copyright (c) 2015 Hantian Zhang, Xiaolin Guo. All rights reserved.


#ifndef LANCZOSINVERSESPARSE
#define LANCZOSINVERSESPARSE

void LanczosInverseSparse(double *A_values, int *A_colIdx, int *A_rowStart, 
    double *M_values, int *M_colIdx, int *M_rowStart, 
    int N,int M_dim,double **V, double **W, 
    double **T_lanc,int qr_iter_num, double *eigenvalues, 
    double **ritz, double threshold);

#endif
