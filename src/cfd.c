#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

// Include grid first to define node,cell, and face types
#include "math_helpers.h"
#include "grid.h"
#include "cfd.h"

/* Indexing macros (row major)*/
//#define IDX(i,j,nx) ((j)*(nx) + (i))
//#define IDX(i,j,eq,nx,ny) ((j)*(nx) + (i) + (eq)*(nx*ny))

#define IDX(i,eq,NCELLS) ( (i) + (eq*NCELLS) )
#define vecIDX(i,j,nx) ( (j)*(nx) + (i) ) //indexing for vectors stored in row major format as [x1,x2,...,xn,y1,y2,...yn,z1,z2,...zn]. nx is number of cells in this case

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

	int NPOINTS = 0, NCELLS = 0, CELL_LIST_SIZE = 0, MAX_FACES = 0, NFACES = 0, NDEGEN_CELLS=0, NBOUNDARIES = 4;


	// Load grid from file and store in nodes and cells arrays, also calculate MAX_FACES for memory allocation of faces array
	int err = read_grid(filename, &nodes, &cells, &NPOINTS, &NCELLS, &CELL_LIST_SIZE, &MAX_FACES, &NDEGEN_CELLS);
	
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

	// Initialize boundaries (change to allocate for more complex gemoetry)
	boundary boundaries[4]; // boundaries
	boundaryType p1_boundaries[4] = { Neumann, Robin, Neumann, Dirichlet };

	for (int i = 0; i < NBOUNDARIES; i++)
	{
		int endpoints[2];
		endpoints[0] = i;
		endpoints[1] = i+1 % (NBOUNDARIES); 
		build_boundary(&boundaries[i], i, endpoints, p1_boundaries[i], nodes, faces, &NFACES);
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

int compute_lsq_gradient(node* nodes, cell* cells, face* faces, int* NCELLS,
	int *NDEGEN_CELLS, int* NFACES, double* phi, double* grad)

{	// Initialize gradient coefficint matrix to zero (See eq. 9.27)
	// Number of cells w/ volume
	int NVOL_CELLS = (*NCELLS) - (*NDEGEN_CELLS);

	// A12 and A21 are the same matrix so only need to allocate one of them
	double* A11 = malloc(NVOL_CELLS * sizeof(double));
	if (A11 == NULL)
	{
		fprintf(stderr, "Error: Memory allocation failed for A11 array.\n");
		return 1;
	}

	double* A12 = malloc(NVOL_CELLS * sizeof(double));
	if (A12 == NULL)
	{
		fprintf(stderr, "Error: Memory allocation failed for A12 array.\n");
		free(A11); // Free previously allocated A11 before exiting
		return 1;
	}

	double* A22 = malloc(NVOL_CELLS * sizeof(double));
	if (A22 == NULL)
	{
		fprintf(stderr, "Error: Memory allocation failed for A22 array.\n");
		free(A11); // Free previously allocated A11 before exiting
		free(A12); // Free previously allocated A12 before exiting
		return 1;
	}

	// Initialize B vector for least squares gradient calculation, size is NCELLS x 2 (x and y components)
	double* b = malloc(NVOL_CELLS * sizeof(double) * 2);
	if (b == NULL)
	{
		fprintf(stderr, "Error: Memory allocation failed for b1 array.\n");
		free(A11); // Free previously allocated A11 before exiting
		free(A12); // Free previously allocated A12 before exiting
		free(A22); // Free previously allocated A22 before exiting
		return 1;
	}


	// Initialize coefficients to zero
	memset(A11, 0, NVOL_CELLS * sizeof(double));
	memset(A12, 0, NVOL_CELLS * sizeof(double));
	memset(A22, 0, NVOL_CELLS * sizeof(double));

	// Loop over all faces and calculate contributions to gradient coefficient matrices
	for (int i = 0; i < *NFACES; i++)
	{
		face* f = &faces[i];

		if (f->neighbor == -1)
		{
			// Boundary face, skip for now (will need to apply boundary conditions later)
			continue;
		}

		cell* C = &cells[f->owner];
		cell* F = &cells[f->neighbor];

		//Define rCF 
		double dxk = f->xc - C->xc;
		double dyk = f->yc - C->yc;
		//double dzk = f->zc - C->zc;

		double dphi = phi[IDX(F->id, 0, *NCELLS)] - phi[IDX(C->id, 0, *NCELLS)];

		// Compute Weight
		double w = 1.0 / sqrt(dxk * dxk + dyk * dyk);

		// Update A11, A12, A22, and b for owner cell
		A11[C->id] += w * dxk * dxk;
		A12[C->id] += w * dxk * dyk; // same as A21
		A22[C->id] += w * dyk * dyk;

		// Update b for owner cell
		b[vecIDX(C->id, 0, NVOL_CELLS)] += w * dphi * dxk; // x component/row of b
		b[vecIDX(C->id, 1, NVOL_CELLS)] += w * dphi * dyk; // y component/row of b

		// Update A11, A12, A22, and b for neighbor cell
		A11[F->id] += w * dxk * dxk;
		A12[F->id] += w * dxk * dyk; // same as A21
		A22[F->id] += w * dyk * dyk;

		// Update b for neighbor cell
		b[vecIDX(F->id, 0, NVOL_CELLS)] += w * dphi * dxk; // x component/row of b
		b[vecIDX(F->id, 1, NVOL_CELLS)] += w * dphi * dyk; // y component/row of b	
	}

	// Loop over all cells and solve for gradient
	for (int i = 0; i < NVOL_CELLS; i++)
	{
		int cell_id = i + *NDEGEN_CELLS - 1; // Adjust index to account for degenerate cells at the beginning of the cells array

		/*
		solve_2x2_system(A11[i], A12[i], A12[i], A22[i], 
			b[vecIDX(i,0,NVOL_CELLS)], b[vecIDX(i, 1, NVOL_CELLS)],
			&grad[])*/
	}

	// Free allocated memory for gradient coefficient matrices
	free(A11);
	free(A12);
	free(A22);
	free(b);
	return 0;
}
