/* ---------------------------------------------
* cfd.h
* Header for CFD solver functions and output
* ---------------------------------------------*/

#ifndef CFD_H
#define CFD_H

#include "grid.h"

// Output VTK file for visualization in ParaView
int write_vtk_output(const char *out_filename, const node *nodes, const cell *cells,
					 const int *NPOINTS, const int *NCELLS, const int *CELL_LIST_SIZE,
					 double **phi, double **grad, const double *p, const double *grad_p);





#endif // !