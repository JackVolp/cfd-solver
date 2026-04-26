#ifndef SOLVER_H
#define SOLVER_H

#include <stdio.h>
#include <math.h>
#include "grid.h"
#include "setup.h"
#include "math_helpers.h"



/* Indexing macros (row major)*/
//#define IDX(i,j,nx) ((j)*(nx) + (i))
//#define IDX(i,j,eq,nx,ny) ((j)*(nx) + (i) + (eq)*(nx*ny))

#define IDX(i,j,nx) ((i) + ((j) * (nx))) // column major for normal arrays
//#define vecIDX(i,j,nx) ( (j)*(nx) + (i) ) //indexing for vectors stored in row major format as [x1,x2,...,xn,y1,y2,...yn,z1,z2,...zn]. nx is number of cells in this case. All the x-dirs, then all the y dirs, then all the z-dirs. I is the index of phi and j is the index of the direction (0,1,2 for x,y,z respectively)

/* -------------------------------------------------------------------------- */
/* Gradient Routines */
/* -------------------------------------------------------------------------- */
// Compute Gradient 
int compute_lsq_gradient(node* nodes, cell* cells, face* faces, int* NCELLS, int* NDEGEN_CELLS, int* NFACES, double* phi, double* grad);

//Interpolate Gradient
int grad2face(double* grad_face, double* grad_C, double* grad_F, double* rCF, double dCF, double phi_C, double phi_F, cell* cell_C, cell* cell_F, face* f);

/* -------------------------------------------------------------------------- */
/* Diffusion Routines */
/* -------------------------------------------------------------------------- */
int build_diffusion(double* A, double* b, double* phi, double* grad, node* nodes, cell* cells, face* faces, boundary* boundaries, int* NCELLS, int* NDEGEN_CELLS, int* NFACES);

/* -------------------------------------------------------------------------- */
/* Advection Routines */
/* -------------------------------------------------------------------------- */
int build_advection(double* A, double* b, double* phi, double* grad, node* nodes, cell* cells, face* faces, boundary* boundaries, int* NCELLS, int* NDEGEN_CELLS, int* NFACES);

double phi2face(double phi_owner,
	double phi_neighbor,
	double mdot_f,
	double* grad,
	cell* owner,
	cell* neighbor,
	face* f,
	advectionScheme scheme);
/* -------------------------------------------------------------------------- */
/* Transient Routines */
/* -------------------------------------------------------------------------- */
int calc_time_step(cell* cells, double* A, int* NCELLS, int* NDEGEN_CELLS, double* t, double* min_dt);

int build_transient(double* A, double* b, double* phi, cell* cells, int* NCELLS, int* NDEGEN_CELLS, double dt);

int explicit_update(double* A, double* b, double* phi, cell* cells, face* faces, double* dt, int* NCELLS, int* NDEGEN_CELLS);

/* -------------------------------------------------------------------------- */
/* Stopping Criteria Routines */
/* -------------------------------------------------------------------------- */
int maxChng(double* phi, double* phi_old, int* NCELLS, int* NDEGEN_CELLS, double* epsilon);

int calc_Residual(double* A,
	double* b,
	double* phi,
	cell* cells,
	face* faces,
	int* NCELLS,
	int* NDEGEN_CELLS,
	int* NFACES,
	double* residual);

/* -------------------------------------------------------------------------- */
/* Boundary routines */
/* -------------------------------------------------------------------------- */
/*
Function to initialize boundary condition
For dirichelet, only phi_b used
for neumann, only q_b used
for robin phi_b and h_infty used
*/
int applyBoundary(boundary* b, cell* cells,
	face* faces, double* phi, double* grad, int* NCELLS);



#endif // !SOLVER_H