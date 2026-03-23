#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

// Include grid first to define node,cell, and face types
#include "grid.h"
#include "cfd.h"

/* Indexing macros (row major)*/
//#define IDX(i,j,nx) ((j)*(nx) + (i))
//#define IDX(i,j,eq,nx,ny) ((j)*(nx) + (i) + (eq)*(nx*ny))

#define IDX(i,eq,NCELLS) ( (i) + (eq*NCELLS) )

int main(void)
{	
	/*----------------------Grid Input Filename------------------*/
	//const char* filename = "C:\\Users\\jtvol\\Documents\\ME696\\Convection-Diffusion\\out\\build\\x64-Debug\\gmsh_grid.vtk"; //Home PC path
	const char* filename = "C:\\Users\\jvolponi0552\\Documents\\GitHub\\cfd-solver\\gmsh_grid.vtk"; //Lab PC path
	/*-----------------------------------------------------------*/


	// Load grid
	node* nodes;
	cell* cells; 
	face* faces;

	int NPOINTS = 0, NCELLS = 0, CELL_LIST_SIZE = 0, MAX_FACES = 0, NFACES;

	// Load grid from file and store in nodes and cells arrays, also calculate MAX_FACES for memory allocation of faces array
	int err = read_grid(filename, &nodes, &cells, &NPOINTS, &NCELLS, &CELL_LIST_SIZE, &MAX_FACES);
	
	if (err != 0)
	{
		fprintf(stderr, "read_grid failed with error code %d\n", err);
		return 1;
	}

	// Calculate Cell Centroid, Volume, Face information, and other geometric properties
	err = build_faces_and_cells(nodes, cells, &NCELLS, &MAX_FACES, &NFACES,&faces);

	int NEQNS = 1; // Number of transport equations solved

	// Allocate conservative scalars
	double* phi = malloc((NEQNS * NCELLS) * sizeof(double));
	if (phi == NULL)
	{
		// Print error message to stderr stream and exit
		fprintf(stderr, "Error: Memory allocation failed for phi array.\n");
		return 1; // Exit with error code
	}

	// initialize phi[eq=0] to 1 everywhere and 0 for other equations
	memset(phi, 0, (NEQNS * NCELLS) * sizeof(double));
	for (int i = 0; i < NCELLS; i++)
	{
		phi[IDX(i, 0, NCELLS)] = 1;
	}

	/* Iteration loop
	* //Calculate gradient at faces
	*	loop over interior faces
	*		find flux at each face
	*		add flux to gradient of owner cell and subtract flux from gradient of neighbor cell
	*	
	*	loop over boundary faces
	*		find flux at each face
	*		add flux to gradient of owner cell
	* 
	*	loop over cells
	*		divide gradient by volume of cell to get final gradient value
	* 
	* //Construct Matrix
	*	initialize aC and aF as zeros for all cells
	*	initialize b as Qc*Vc for all cells
	*	
	*	loop over all internal faces
	*		af(owner) +=  -gamma(face)*Ef(face)/dCF
	*		ac(owner) +=  gamma(face)*Ef(face)/dCF
	* 
	*		af(neighbor) +=  gamma(face)*Ef(face)/dCF
	*		ac(neighbor) +=  -gamma(face)*Ef(face)/dCF
	* 
	*		interpolate gradient to face
	*		
	*		b(owner) += gamma(face)*grad(face) dot Tf(face)
	*		b(neighbor) += gamma(face)*grad(face) dot -Tf(face)
	* 
	*	loop over all boundary faces
	* 
	*		switch (boundary condition type)
	*		case Dirichlet: 
	*			dirichelet coefficents pg 252
	*		case Neumann:
	*			neumann coefficients pg 221
	*		case Robin:
	*			robin/mixed coefficients pg 253-254

	// Apply Boundary conditions

	// Calculate gradient 
	// loop over interior faces
	//	find flux at faces 
	//	add flux to owner cell and subtract flux from neighbor cell
	//	loop over boundary faces
	//  calculate boundary flux at faces and add to owner cell 
	//	loop over all elemetns and divide gardient by voume of element
	// 
	//  
	*/

	// ------ Write output file --------
	err = write_vtk_output("output_file.vtk", &nodes, &cells, &NPOINTS, &NCELLS,
		&CELL_LIST_SIZE, &phi);
	if (err != 0)
	{
		fprintf(stderr, "write_vtk_output failed with error code %d\n", err);
		return 1;
	}

	// Release Allocated Memory for grid
	free_grid(nodes, cells, faces, NCELLS,NFACES);

	// Release conservative scalars memory
	free(phi);

	

	printf("To C or not to C: that is the question. \n");
	return 0;
	
}


/*---------------------------------------------------------------------------
* Write data function and grid 
----------------------------------------------------------------------------*/
int write_vtk_output(const char* out_filename,	node** nodes,	cell** cells,
	int* NPOINTS,	int* NCELLS,	int* CELL_LIST_SIZE,	double** phi)
{
	//Open file for writing
	FILE* fp = fopen(out_filename, "w");

	if (!fp)
	{
		perror("Error writing data file\n");
		return 1;
	}

	// Header section
	fprintf(fp, "# vtk DataFile Version 3.0\n");
	fprintf(fp, "Jack's first unstructured vtk datafile\n");
	fprintf(fp, "ASCII\n");

	// Dataset type definition
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

	// Write points
	fprintf(fp, "POINTS %d double\n", *NPOINTS);
	
	for (int i = 0; i < *NPOINTS; i++) 
	{
		fprintf(fp, "%.15f %.15f %.15f\n", (*nodes)[i].x, (*nodes)[i].y, (*nodes)[i].z);
	}

	// Write Cells
	fprintf(fp, "CELLS %d %d\n", *NCELLS, *CELL_LIST_SIZE);
	
	for (int i = 0; i < *NCELLS; i++)
	{
		// Get number of nodes for current cell
		int num_nodes = (*cells)[i].num_nodes;
		fprintf(fp, "%d ", num_nodes);

		for (int point = 0; point < num_nodes; point++)
		{
			if (point == (num_nodes - 1))
			{
				fprintf(fp, "%d\n", (*cells)[i].node_ids[point]);
			}
			else
			{
				fprintf(fp, "%d ", (*cells)[i].node_ids[point]);
			}
			
		}
	}

	// Write Cell Types
	fprintf(fp, "CELL_TYPES %d \n", *NCELLS);

	for (int i = 0; i < *NCELLS; i++)
	{
		fprintf(fp, "%d\n", (*cells)[i].type);
	}

	// Write Scalar data 
	fprintf(fp, "CELL_DATA %d\n", *NCELLS);

	// Write cell id
	fprintf(fp, "SCALARS cellId int 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");

	for (int i = 0; i < *NCELLS; i++)
	{
		fprintf(fp, "%d\n", (*cells)[i].id);
	}
	// Write phi
	fprintf(fp, "SCALARS phi[%d] double 1\n",0);
	fprintf(fp, "LOOKUP_TABLE default\n");

	for (int i = 0; i < *NCELLS; i++)
	{
		fprintf(fp, "%.15f\n", (*phi)[IDX(i, 0, *NCELLS)]);
	}

	// Write Node data
	fprintf(fp, "POINT_DATA %d\n", *NPOINTS);
	// Write node id
	fprintf(fp, "SCALARS nodeId int 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");

	for (int i = 0; i < *NPOINTS; i++)
	{
		fprintf(fp, "%d\n", (*nodes)[i].id);
	}

	fclose(fp);

	return 0;
}

int compute_gradient(node* nodes, cell* cells, face* faces, int* NCELLS,
	int* NFACES, double* phi, double* grad)
{
	return 0;
}
