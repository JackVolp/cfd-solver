#ifndef SOLVER_H
#define SOLVER_H

#include <stdio.h>
#include <math.h>
#include "grid.h"
#include "constants.h"
#include "math_helpers.h"

/* Indexing macros (row major)*/
//#define IDX(i,j,nx) ((j)*(nx) + (i))
//#define IDX(i,j,eq,nx,ny) ((j)*(nx) + (i) + (eq)*(nx*ny))

#define IDX(i,j,nx) ( (i) + (j*nx) ) // column major for normal arrays
//#define vecIDX(i,j,nx) ( (j)*(nx) + (i) ) //indexing for vectors stored in row major format as [x1,x2,...,xn,y1,y2,...yn,z1,z2,...zn]. nx is number of cells in this case. All the x-dirs, then all the y dirs, then all the z-dirs. I is the index of phi and j is the index of the direction (0,1,2 for x,y,z respectively)


int compute_lsq_gradient(node* nodes, cell* cells, face* faces, int* NCELLS, int* NDEGEN_CELLS, int* NFACES, double* phi, double* grad);

int build_matrix(double* A, double* b, double* phi, double* grad, node* nodes, cell* cells, face* faces, boundary* boundaries, int* NCELLS, int* NDEGEN_CELLS, int* NFACES);

int grad2face(double* grad_face, double* grad_C, double* grad_F, double* rCF, double dCF, double phi_C, double phi_F,cell* cell_C, cell* cell_F, face* f);

int maxChng(double* phi, double* phi_old, int* NCELLS, int* NDEGEN_CELLS, double* epsilon);

/*
Function to initialize boundary condition
For dirichelet, only phi_b used
for neumann, only q_b used
for robin phi_b and h_infty used
*/
int applyBoundary(boundary* b, cell* cells,
	face* faces, double* phi, double* grad, int* NCELLS);

#endif // !SOLVER_H