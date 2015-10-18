//  "Generalized-restart Lanczos Method for Schrödinger Equation in Semiconductor Simulation"
//  Created by Hantian Zhang, Xiaolin Guo on 10/09/15.
//  Copyright (c) 2015 Hantian Zhang, Xiaolin Guo. All rights reserved.

#ifndef __QR__
#define __QR__

void qr(double ** A, int M, int N,double * eigenvalues, double **eigenvectors,int niters);

#endif
