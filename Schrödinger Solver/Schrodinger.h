//  "Generalized-restart Lanczos Method for Schrödinger Equation in Semiconductor Simulation"
//  "Section 3: Extraction of m lowest eigenstates"
//  Created by Hantian Zhang, Xiaolin Guo on 10/09/15.
//  Copyright (c) 2015 Hantian Zhang, Xiaolin Guo. All rights reserved.

#ifndef SCHRODINGER
#define SCHRODINGER

void Schrodinger(int N_Sch, int Nx_Sch, int Ny_Sch,
                      int M_dim, int qr_iter,
                      double THRESHOLD, double THRESHOLD_inverse,
                      double *Domain);

#endif