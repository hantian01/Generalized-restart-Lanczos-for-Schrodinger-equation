Schrödinger equation solver in finite element formalism by generalized-restart Lanczos method
===
* A MPI-paralleled solver for 2D time-independent Schrödinger equation in paper: "Generalized-restart Lanczos Method for Schrödinger Equation in Semiconductor Simulation" by Hantian Zhang and Xiaolin Guo.
* Numerical test of "Section 4.2.1: Test on 2D quantum harmonic oscillator".
* Note the potential is modified according to "Section 4.3.1: Scaling".

Feature
--------
* Fast compute several quantum bounded states simultaneously of large-scale semiconductor simulation by Finite Element Method
* Generalized-restart Lanczos iteration and inverse variant employed
* CG_mpi as parallel linear system solver


Usage
--------
compile with the default Makefile, and run with

	mpirun -n procs ./test_mpi Nx

where we take Nx = Ny and procs is the number of MPI processors. Please make sure (Nx^2) % procs == 0.

