//  "Generalized-restart Lanczos Method for Schrödinger Equation in Semiconductor Simulation"
//  Created by Hantian Zhang, Xiaolin Guo on 10/09/15.
//  Copyright (c) 2015 Hantian Zhang, Xiaolin Guo. All rights reserved.

#include "LanczosSparse.h"
#include "LanczosInverseSparse.h"
#include "Sch_Matrix.h"
#include "util.h"
#include "Schrodinger.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

#include <mpi.h>


using std::cout;
using std::endl;

#ifdef TIMER
#include "timer.h"
double timer_CG;
double timer_CG_mpi;
double timer_smvm;
double timer_septdiagmvm;
double timer_ghostmvm;
int counter;
int GS_counter;
int CG_counter;
int InvCounter;
#endif


int main(int argc, char *argv[])
{
#ifdef TIMER
    timer_smvm = 0.;
    timer_CG_mpi = 0.;
#endif
    double timer_total =  0.0;

    int Nx_Sch = 80;
    if (argc == 2)
        Nx_Sch = atoi(argv[1]); // X-axis mesh nodes input
    int Ny_Sch = Nx_Sch;
    int N_Sch = Nx_Sch * Ny_Sch; // interior mesh nodes for Solving Schrödinger equation

    int M_dim = 5; // Subspace dimension, m-th lowest eigenpairs to be computed
    double THRESHOLD = 0.001; // threshold for normal lanczos
    double THRESHOLD_inverse = 0.000001; // threshold for inverse lanczos
    int qr_iter = 60; // QR iteration cycles

    // Domain: [0]->a, [1]->b: Schrodinger Domain
    double Domain[4];
    Domain[0] = -10.; Domain[1] = 10.;
    Domain[2] = -10.; Domain[3] = 10.;

    /*            ||
    || MPI starts ||
    ||            */

    MPI::Init(argc, argv);
    double start = MPI_Wtime();

    int rank, procs;
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        cout << "---- I am running on " << procs << " procs!" << endl;
        cout << "---- Nx = " << Nx_Sch << "  N_Sch = " << N_Sch << endl;

    }

    // Solve Schrödinger equation
    Schrodinger(N_Sch, Nx_Sch, Ny_Sch, M_dim, qr_iter, THRESHOLD, THRESHOLD_inverse, Domain);

    timer_total += MPI_Wtime() - start;

#ifdef TIMER
    if (rank == 0)
    {
        cout << "time for smvm: " << timer_smvm << endl;
        cout << "time for CG_mpi: " << timer_CG_mpi << endl;
        cout << "total time: " << timer_total << endl;
    }
#endif
    

    MPI::Finalize();

    return 0;
}