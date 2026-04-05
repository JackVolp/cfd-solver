#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

// Include grid first to define node,cell, and face types
#include "math_helpers.h"
#include "grid.h"
#include "cfd.h"
#include "constants.h"
#include "math_helpers.h"

/* Indexing macros (row major)*/
//#define IDX(i,j,nx) ((j)*(nx) + (i))
//#define IDX(i,j,eq,nx,ny) ((j)*(nx) + (i) + (eq)*(nx*ny))

#define IDX(i,eq,NCELLS) ( (i) + (eq*NCELLS) )
#define vecIDX(i,j,nx) ( (j)*(nx) + (i) ) //indexing for vectors stored in row major format as [x1,x2,...,xn,y1,y2,...yn,z1,z2,...zn]. nx is number of cells in this case. All the x-dirs, then all the y dirs, then all the z-dirs. I is the index of phi and j is the index of the direction (0,1,2 for x,y,z respectively)

int main(void)
{	
	/*----------------------Grid Input Filename------------------*/
	const char* filename = "C:\\Users\\jtvol\\Documents\\ME696\\Convection-Diffusion\\out\\build\\x64-Debug\\gmsh_grid.vtk"; //Home PC path
	//const char* filename = "C:\\Users\\jvolponi0552\\Documents\\GitHub\\cfd-solver\\gmsh_grid.vtk"; //Lab PC path
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

	// Allocate gradient array (3 components for x,y,z)
	double* grad = malloc((3 * NCELLS) * sizeof(double));
	if (grad == NULL)
	{
		// Print error message to stderr stream and exit
		fprintf(stderr, "Error: Memory allocation failed for grad array.\n");
		return 1; // Exit with error code
	}

	// initialize phi[eq=0] to 1 everywhere and 0 for other equations
	memset(phi, 0, (NEQNS * NCELLS) * sizeof(double));
	for (int i = 0; i < NCELLS; i++)
	{
		phi[IDX(i, 0, NCELLS)] = 1;
	}

	// initialize grad to zero
	memset(grad, 0, ((int)3 * NCELLS) * sizeof(double));

	// initialize gradient array
	// Initialize boundaries (change to allocate for more complex gemoetry)
	boundary boundaries[4]; // boundaries
	boundaryType p1_boundaries[4] = { Neumann, Robin, Neumann, Dirichlet };
	boundaryData p1_boundary_data[4] = {
		{.q_b = 0.0},
		{.robin = {.h_inf = 100.0, .phi_inf = 25.}},
		{.q_b = 0.0},
		{.phi_b = 100.}
	};

	for (int i = 0; i < NBOUNDARIES; i++)
	{
		int endpoints[2];
		endpoints[0] = i;
		endpoints[1] = (i+1) % (NBOUNDARIES); 
		build_boundary(&boundaries[i], i, endpoints, p1_boundaries[i],p1_boundary_data[i], nodes, faces, &NFACES);
	}
	
	// Initialize phi on the boundaries
	for (int i = 0; i < NBOUNDARIES; i++)
	{
		err = initBoundary(&boundaries[i], cells, faces, phi, grad, &NCELLS);

		if (err != 0)
		{
			fprintf(stderr, "initBoundary failed with error code %d\n", err);
			return 1;
		}
	}

	// Solver Loop
	for (int i = 0; i < MAX_ITER; i++)
	{
		err = compute_lsq_gradient(nodes, cells, faces, &NCELLS, &NDEGEN_CELLS, &NFACES, phi, grad);
		if (err != 0)
		{
			fprintf(stderr, "compute_lsq_gradient failed with error code %d\n", err);
			return 1;
		}
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
		&CELL_LIST_SIZE, &phi, &grad);
	if (err != 0)
	{
		fprintf(stderr, "write_vtk_output failed with error code %d\n", err);
		return 1;
	}

	// Release Allocated Memory for grid
	free_grid(nodes, cells, faces, NCELLS,NFACES);

	// Release conservatiWve scalars memory
	free(phi);
	free(grad);

	

	printf("To C or not to C: that is the question. \n");
	return 0;
	
}


/*---------------------------------------------------------------------------
* Write data function and grid 
----------------------------------------------------------------------------*/
int write_vtk_output(const char* out_filename, node** nodes, cell** cells,
	int* NPOINTS, int* NCELLS, int* CELL_LIST_SIZE, double** phi, double** grad)
{
	// Open file for writing
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
		fprintf(fp, "%.15f %.15f %.15f\n",
			(*nodes)[i].x, (*nodes)[i].y, (*nodes)[i].z);
	}

	// Write cells
	fprintf(fp, "CELLS %d %d\n", *NCELLS, *CELL_LIST_SIZE);
	for (int i = 0; i < *NCELLS; i++)
	{
		int num_nodes = (*cells)[i].num_nodes;
		fprintf(fp, "%d ", num_nodes);

		for (int point = 0; point < num_nodes; point++)
		{
			if (point == num_nodes - 1)
			{
				fprintf(fp, "%d\n", (*cells)[i].node_ids[point]);
			}
			else
			{
				fprintf(fp, "%d ", (*cells)[i].node_ids[point]);
			}
		}
	}

	// Write cell types
	fprintf(fp, "CELL_TYPES %d\n", *NCELLS);
	for (int i = 0; i < *NCELLS; i++)
	{
		fprintf(fp, "%d\n", (*cells)[i].type);
	}

	// Write cell data for ALL cells, including degenerate cells
	fprintf(fp, "CELL_DATA %d\n", *NCELLS);

	// Cell IDs
	fprintf(fp, "SCALARS cellId int 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (int i = 0; i < *NCELLS; i++)
	{
		fprintf(fp, "%d\n", (*cells)[i].id);
	}

	// Phi
	fprintf(fp, "SCALARS phi[%d] double 1\n", 0);
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (int i = 0; i < *NCELLS; i++)
	{
		fprintf(fp, "%.15f\n", (*phi)[IDX(i, 0, *NCELLS)]);
	}

	// Gradient vector
	fprintf(fp, "VECTORS grad_phi double\n");
	for (int i = 0; i < *NCELLS; i++)
	{
		fprintf(fp, "%.15f %.15f %.15f\n",
			(*grad)[vecIDX(i, 0, *NCELLS)],
			(*grad)[vecIDX(i, 1, *NCELLS)],
			(*grad)[vecIDX(i, 2, *NCELLS)]);
	}

	// Optional: mark degenerate cells explicitly
	fprintf(fp, "SCALARS isDegenerate int 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (int i = 0; i < *NCELLS; i++)
	{
		fprintf(fp, "%d\n", ((*cells)[i].volume <= 0.0) ? 1 : 0);
	}

	// Point data
	fprintf(fp, "POINT_DATA %d\n", *NPOINTS);
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
	memset(b, 0, NVOL_CELLS * sizeof(double) * 2);

	// Loop over all faces and calculate contributions to gradient coefficient matrices
	for (int i = 0; i < *NFACES; i++)
	{
		face* f = &faces[i];

		//if (f->neighbor == -1)
		//{
		//	// Boundary face, skip for now (will need to apply boundary conditions later)
		//	continue;
		//}

		cell* C = &cells[f->owner];
		cell* F = &cells[f->neighbor];

		// Indicies for arrays with size of volume cells 
		int vC_idx = C->id - *NDEGEN_CELLS;
		int vF_idx = F->id - *NDEGEN_CELLS;

		//Define rCF 
		double dxk = f->xc - C->xc;
		double dyk = f->yc - C->yc;
		//double dzk = f->zc - C->zc;

		double dphi = phi[IDX(F->id, 0, *NCELLS)] - phi[IDX(C->id, 0, *NCELLS)];

		// Compute Weight
		double w = 1.0 / sqrt(dxk * dxk + dyk * dyk);

		// Update A11, A12, A22, and b for owner cell
		A11[vC_idx] += w * dxk * dxk;
		A12[vC_idx] += w * dxk * dyk; // same as A21
		A22[vC_idx] += w * dyk * dyk;

		// Update b for owner cell
		b[vecIDX(vC_idx, 0, NVOL_CELLS)] += w * dphi * dxk; // x component/row of b
		b[vecIDX(vC_idx, 1, NVOL_CELLS)] += w * dphi * dyk; // y component/row of b

		// Update A11, A12, A22, and b for neighbor cell if the neighbor cell is not degenerate. Neighbor will be degenerate for boundary cells.
		if (!f->boundary_face)
		{
			A11[vF_idx] += w * dxk * dxk;
			A12[vF_idx] += w * dxk * dyk; // same as A21
			A22[vF_idx] += w * dyk * dyk;

			// Update b for neighbor cell
			b[vecIDX(vF_idx, 0, NVOL_CELLS)] += w * dphi * dxk; // x component/row of b
			b[vecIDX(vF_idx, 1, NVOL_CELLS)] += w * dphi * dyk; // y component/row of b
		}
		

			
	}

	// Loop over all cells and solve for gradient
	for (int i = 0; i < NVOL_CELLS; i++)
	{
		int cell_id = i + *NDEGEN_CELLS; // Adjust index to account for degenerate cells at the beginning of the cells array

		solve_2x2_system(A11[i], A12[i], A12[i], A22[i], 
			b[vecIDX(i,0,NVOL_CELLS)], b[vecIDX(i, 1, NVOL_CELLS)],
			&grad[vecIDX(cell_id, 0, *NCELLS)], &grad[vecIDX(cell_id, 1, *NCELLS)]); // Store gradient in correct location in grad array based on cell id
	}

	// Free allocated memory for gradient coefficient matrices
	free(A11);
	free(A12);
	free(A22);
	free(b);
	return 0;
}

int initBoundary(boundary* b, cell* cells, 
	face* faces, double* phi, double* grad, int* NCELLS)
{
	// Loop over all faces in boundary and apply boundary conditions
	for (int i = 0; i < b->num_faces; i++)
	{
		int face_id = b->face_ids[i];
		face* f = &faces[face_id];
		cell* c_owner = &cells[f->owner];

		// Should always be a boundary face but just to be sure
		int phi_face_idx = f->neighbor ; // If boundary face, use owner cell for phi index, otherwise use neighbor cell (should not trigger for internal faces)

		int phi_owner_idx = f->owner; // Index for owner cell in phi array

		switch (b->type)
		{
		case Dirichlet: 
			// For Dirichlet, we can set the boundary value directly
			phi[phi_face_idx] = b->data.phi_b; // Set phi at owner cell to boundary value
			break;
		case Neumann: {
			// magnitude of surface area vector
			double Ef[3] = { f->Ex, f->Ey, f->Ez }; // Surface area vector of face
			double Ef_mag = 0;
			magnitude(Ef, &Ef_mag);

			// Distance from owner cell centroid to face centroid
			// Always positive since face vector always points outward from owner cell
			double rCF[3] = { f->xc - c_owner->xc, f->yc - c_owner->yc, f->zc - c_owner->zc }; // vector from cell centroid to face centroid

			double d_CF = 0;
			magnitude(rCF, &d_CF);

			double gDiff = Ef_mag / d_CF; // "Geometric Diffusion Coefficient"

			phi[phi_face_idx] = (GAMMA * gDiff * phi[phi_owner_idx] - b->data.q_b) 
				/ (GAMMA * gDiff); // Eq 8.42 from textbook
			break;
		}
		case Robin: {
			// For Robin, we need to calculate the equivalent boundary value based on the given phi_b, q_b, and h_infty
			// Ignore cross diffusion term for initialization? maybe. Dont have gradient at face yet
			
			double h_inf = b->data.robin.h_inf;
			double phi_inf = b->data.robin.phi_inf;

			// magnitude of surface area vector
			double Ef[3] = { f->Ex, f->Ey, f->Ez }; // Surface area vector of face
			double Ef_mag = 0;
			magnitude(Ef, &Ef_mag);

			// Distance from owner cell centroid to face centroid
			// Always positive since face vector always points outward from owner cell
			double rCF[3] = { f->xc - c_owner->xc, f->yc - c_owner->yc, f->zc - c_owner->zc }; // vector from cell centroid to face centroid

			double d_CF = 0;
			magnitude(rCF, &d_CF);

			double gDiff = Ef_mag / d_CF; // "Geometric Diffusion Coefficient"

			double Sf[3] = { f->Sx, f->Sy, f->Sz }; // Surface area vector of face
			double Sf_mag = 0;
			magnitude(Sf, &Sf_mag);

			// Compute phi at the boundary using eq.8.85	
			double gradPhi_face[3] = { grad[vecIDX(phi_face_idx, 0, *NCELLS)], grad[vecIDX(phi_face_idx, 1, *NCELLS)], grad[vecIDX(phi_face_idx, 2, *NCELLS)] }; 

			double grad_dot_Tf = gradPhi_face[0] * f->Tx + gradPhi_face[1] * f->Ty + gradPhi_face[2] * f->Tz; // Gradient at face dot tangential contribution vector of face

			phi[phi_face_idx] = (h_inf * Sf_mag * phi_inf
				+ GAMMA*gDiff*phi[phi_owner_idx]
				- GAMMA*grad_dot_Tf) / (h_inf*Sf_mag + GAMMA*gDiff);
			break;
		}
		default:
			fprintf(stderr, "Error: Unknown boundary condition type for boundary ID %d\n", b->id);
			return 1;
		}
	}

	
	return 0;
}