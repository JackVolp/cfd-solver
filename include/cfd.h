/* ---------------------------------------------
* cfd.h
* Header for CFD solver functions and output
* ---------------------------------------------*/

#ifndef CFD_H
#define CFD_H

#include "grid.h"

// Output VTK file for visualization in ParaView
int write_vtk_output(const char* out_filename,
	node** nodes,
	cell** cells,
	int* NPOINTS,
	int* NCELLS,
	int* CELL_LIST_SIZE,
	double** scalars);

/* 
Function to initialize boundary condition
For dirichelet, only phi_b used
for neumann, only q_b used
for robin phi_b and h_infty used
*/
int initBoundary(boundary* b, double phi_b, double q_b, double h_infty, cell* cells, face* faces, double* phi);

int compute_lsq_gradient(node* nodes, cell* cells, face* faces, int* NCELLS,
	int* NDEGEN_CELLS, int* NFACES, double* phi, double* grad);

#endif // !