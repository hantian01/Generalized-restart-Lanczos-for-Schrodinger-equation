//  "Generalized-restart Lanczos Method for Schrödinger Equation in Semiconductor Simulation"
//  "Section 4.3: Distributed-memory parallelism"
//  Created by Hantian Zhang, Xiaolin Guo on 10/09/15.
//  Copyright (c) 2015 Hantian Zhang, Xiaolin Guo. All rights reserved.

#ifndef __CG_MPI__
#define __CG_MPI__

#include "util.h"

void CG_mpi_wrapper(int                       N,
                    CSRLocalMatrix_t          &csr_loc,
                    double                    *b,
                    double                    *x,
                    double                    threshold);

void get_local_matrix(double                    *values,
                      int                       *colIdx,
                      int                       *rowStart,
                      int                       N,
                      CSRLocalMatrix_t          &csr_loc);

#endif
