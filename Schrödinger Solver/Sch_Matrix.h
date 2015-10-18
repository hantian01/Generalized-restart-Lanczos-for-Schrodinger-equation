//  "Generalized-restart Lanczos Method for Schrödinger Equation in Semiconductor Simulation"
//  "Section 4.1: Linear FEM on 2D triangular mesh"
//  Created by Hantian Zhang, Xiaolin Guo on 10/09/15.
//  Copyright (c) 2015 Hantian Zhang, Xiaolin Guo. All rights reserved.

#ifndef SCH_MATRIX
#define SCH_MATRIX

// static double hbar = 1.05457172647e-16; // Unit: kg nm^2/s
// static double q = 1.60217656535e-19; // Unit: C
// static double k_b = 1.380648813e-5; // Unit: kg nm^2 / s^2 K
// static double epsilon_0 = 8.8541878176e-12; // Unit: kg / s^2
// static double Tempature = 293.15;

// static double hbar = 1.0; // Unit: kg nm^2/s
// static double q = 1.0; // Unit: C
// static double k_b = 1.0; // Unit: kg nm^2 / s^2 K
// static double epsilon_0 = 1.0; // Unit: kg / s^2
// static double Tempature = 1;

void GalerkinMatrix_Sch(double*& A_values, int*& A_colIdx, int*& A_rowStart,
                        double*& M_values, int*& M_colIdx, int*& M_rowStart,
                        int Nx_Sch, int Ny_Sch,
                        double x_a, double x_b, double y_a, double y_b
                       );

#endif
