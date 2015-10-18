//  "Generalized-restart Lanczos Method for Schrödinger Equation in Semiconductor Simulation"
//  "Section 4.1: Linear FEM on 2D triangular mesh"
//  Created by Hantian Zhang, Xiaolin Guo on 10/09/15.
//  Copyright (c) 2015 Hantian Zhang, Xiaolin Guo. All rights reserved.

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <list>
#include <vector>
#include <algorithm>

#include "Sch_Matrix.h"
#include "util.h"

using namespace std;

// "Section 4.2.1: Test on 2D quantum harmonic oscillator"
inline double V(double x, double y, double m)
{
	double omega = 1.;
	// "Section 4.3.1: Scaling"
	// to avoid degeneracy of eigenvalue, change potential term to the following:
	return 0.5 * m * omega * omega * (x * x + 8 * y * y);
}

inline void indexMap(int Nx, int Ny, int idxEle, int *idxGlo, double vertice[3][2],
                     double hx, double hy, double xa, double ya)
{
	int nx(0), ny(0);
	if (idxEle % 2 != 0)   // odd element index
	{
		ny = idxEle / (2 * (Nx - 1));
		nx = (int) (( idxEle + 1 - 2 * ny * (Nx - 1) ) / 2);
		// if (idxLoc == 0)
		idxGlo[0] = (idxEle + 1) / 2 + ny - 1;
		vertice[0][0] = xa + ( nx - 1 ) * hx; // x position in global coordinate
		vertice[0][1] = ya + ny * hy; // y position in global coordinate
		// else if (idxLoc == 1)
		idxGlo[1] = idxGlo[0] + 1;
		vertice[1][0] = vertice[0][0] + hx; // x position
		vertice[1][1] = vertice[0][1]; // y position
		// else if (idxLoc == 2)
		idxGlo[2] = idxGlo[0] + Nx ;
		vertice[2][0] = vertice[0][0]; // x position
		vertice[2][1] = vertice[0][1] + hy; // y position;
	}
	else  // even element index
	{
		ny = (idxEle - 1) / (2 * (Nx - 1));
		nx = (int) (( idxEle - 2 * ny * (Nx - 1) ) / 2 );
		// if (idxLoc == 0)
		idxGlo[0] = (ny + 1) * Nx + nx + 1 - 1;
		vertice[0][0] = xa + nx * hx; // x position in global coordinate
		vertice[0][1] = ya + ( ny + 1 ) * hy; // y position in global coordinate
		// else if (idxLoc == 1)
		idxGlo[1] = idxGlo[0] - 1;
		vertice[1][0] = vertice[0][0] - hx; // x position
		vertice[1][1] = vertice[0][1]; // y position

		// else if (idxLoc == 2)
		idxGlo[2] = idxGlo[0] - Nx;
		vertice[2][0] = vertice[0][0]; // x position
		vertice[2][1] = vertice[0][1] - hy; // y position
	}
}

void GalerkinMatrix_Sch(double*& A_values, int*& A_colIdx, int*& A_rowStart,
                        double*& M_values, int*& M_colIdx, int*& M_rowStart,
                        int Nx_Sch, int Ny_Sch,
                        double x_a, double x_b, double y_a, double y_b
                       )
{
	int N, N_element;
	double area, hx, hy;
	double m = 1.0;
	double hbar = 1.0;
	double sigma = hbar * hbar / (2 * m);

	int idxGlo[3];
	double a1_loc[3][3], a2_loc[3][3], m_loc[3][3];
	double vertice[3][2];

	N = Nx_Sch * Ny_Sch;
	N_element =  2 * (Nx_Sch - 1) * (Ny_Sch - 1);
	hx = (x_b - x_a) / (Nx_Sch - 1);
	hy = (y_b - y_a) / (Ny_Sch - 1);
	area = hx * hy / 2;

	// Local matrix computation for a1
	double coef1 = sigma / 2 * hx * hy;
	a1_loc[0][0] = coef1 * (1 / (hx * hx) + 1 / (hy * hy));
	a1_loc[0][1] = - coef1 / (hx * hx); a1_loc[1][0] = a1_loc[0][1];
	a1_loc[0][2] = - coef1 / (hy * hy); a1_loc[2][0] = a1_loc[0][2];
	a1_loc[1][1] = coef1 / (hx * hx);
	a1_loc[1][2] = 0; a1_loc[2][1] = 0;
	a1_loc[2][2] = coef1 / (hy * hy);

	// local matrix computation for m
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
			m_loc[i][j] = area / 12;
		m_loc[i][i] *= 2;
	}

	std::vector< std::list<maRowEle> > ma(N), mm(N);
	// global assemble and local computation for a2
	for (int idxEle = 1; idxEle <= N_element; ++idxEle)
	{
		indexMap(Nx_Sch, Ny_Sch, idxEle, idxGlo, vertice, hx, hy, x_a, y_a);

		double V12 = V((vertice[0][0] + vertice[1][0]) / 2, (vertice[0][1] + vertice[1][1]) / 2, m);
		double V13 = V((vertice[0][0] + vertice[2][0]) / 2, (vertice[0][1] + vertice[2][1]) / 2, m);
		double V23 = V((vertice[2][0] + vertice[1][0]) / 2, (vertice[2][1] + vertice[1][1]) / 2, m);

		double coef2 = area / 12;
		a2_loc[0][0] = coef2 * (V13 + V12);
		a2_loc[0][1] = coef2 * V12; a2_loc[1][0] = a2_loc[0][1];
		a2_loc[0][2] = coef2 * V13; a2_loc[2][0] = a2_loc[0][2];
		a2_loc[1][1] = coef2 * (V12 + V23);
		a2_loc[1][2] = coef2 * V23; a2_loc[2][1] = a2_loc[1][2];
		a2_loc[2][2] = coef2 * (V13 + V23);

		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				addToMA(a1_loc[i][j] + a2_loc[i][j], idxGlo[i], idxGlo[j], ma);
				addToMA(m_loc[i][j], idxGlo[i], idxGlo[j], mm);
			}
		}
	}

	// Convert to CSR form
	maToCSR(ma, A_values, A_colIdx, A_rowStart);
	maToCSR(mm, M_values, M_colIdx, M_rowStart);


	// printCSR(N, M_values, M_colIdx, M_rowStart);
	// printCSR(N, A_values, A_colIdx, A_rowStart);

}
