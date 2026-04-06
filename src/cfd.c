#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

// Include grid first to define node,cell, and face types
#include "math_helpers.h"
#include "grid.h"
#include "cfd.h"
#include "solver.h"
#include "constants.h"
#include "math_helpers.h"


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

	int NPOINTS = 0, NCELLS = 0, CELL_LIST_SIZE = 0, MAX_FACES = 0, NFACES = 0, NDEGEN_CELLS=0, NBOUNDARIES = 4, NSOLCELLS = 0;


	// Load grid from file and store in nodes and cells arrays, also calculate MAX_FACES for memory allocation of faces array
	int err = read_grid(filename, &nodes, &cells, &NPOINTS, &NCELLS, &CELL_LIST_SIZE, &MAX_FACES, &NDEGEN_CELLS);

	NSOLCELLS = NCELLS - NDEGEN_CELLS; // Number of cells that have volume and are included in the solution
	
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

	lapack_int n = NSOLCELLS; // Number of equations (size of the system)
	lapack_int nrhs = 1; // Number of right-hand sides (columns of B)
	lapack_int lda = NSOLCELLS; // Leading dimension of A
	lapack_int ldb = NSOLCELLS; // Leading dimension of B

	lapack_int* ipiv = malloc(NSOLCELLS * sizeof(lapack_int)); // Pivot indices for LU factorization
	if (!ipiv)
	{
		// Print error message to stderr stream and exit
		fprintf(stderr, "Error: Memory allocation failed for ipiv array.\n");
		return 1; // Exit with error code
	}

	//
	// Allocate and initialize matrix coefficients and source term vector for linear system
	// For each row (cell):
	// one aC item
	// num_faces aF items
	// replace with sparse storage later 
	//
	//double* aC = malloc(NCELLS * sizeof(double)); // Diagonal coefficients
	//if (!aC)
	//{
	//	// print error message to stderr stream
	//	fprintf(stderr, "Error: Memory allocation failed for aC array.\n");
	//	return 1; // Exit with error code
	//}

	//double* aF = malloc((MAX_FACES) * sizeof(double)); // Off-diagonal coefficients (face contributions)
	//if (!aF)
	//{
	//	// print error message to stderr stream
	//	fprintf(stderr, "Error: Memory allocation failed for aF array.\n");
	//	return 1; // Exit with error code
	//}

	//int* f_col = malloc((MAX_FACES) * sizeof(int)); // Column indices for off-diagonal coefficients (cell id of the neighbor cells
	//if (!f_col)
	//{
	//	// print error message to stderr stream
	//	fprintf(stderr, "Error: Memory allocation failed for f_col array.\n");
	//	return 1; // Exit with error code
	//}
	double* A = malloc((NEQNS * NSOLCELLS * NSOLCELLS) * sizeof(double)); // Coefficient matrix (will be stored in sparse format later)
	if (!A)
	{
		// print error message to stderr stream
		fprintf(stderr, "Error: Memory allocation failed for A array.\n");
		return 1; // Exit with error code
	}

	double* b = malloc(NSOLCELLS * sizeof(double)); // Source term vector
	if (!b)
	{
		// print error message to stderr stream
		fprintf(stderr, "Error: Memory allocation failed for b array.\n");
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

	// ------------Solver Loop---------------------------
	for (int i = 0; i < MAX_ITER; i++)
	{
		// Compute Gradient
		err = compute_lsq_gradient(nodes, cells, faces, &NCELLS, &NDEGEN_CELLS, &NFACES, phi, grad);
		if (err != 0)
		{
			fprintf(stderr, "compute_lsq_gradient failed with error code %d\n", err);
			return 1;
		}

		// Build Matrix and Source term
			// initialize matrix coefficients and source term vector to zero
		memset(A, 0, (NEQNS* NSOLCELLS* NSOLCELLS) * sizeof(double));
		memset(b, 0, ((NEQNS* NSOLCELLS) * sizeof(double)));
		err = build_matrix(A, b, phi, grad, nodes, cells, faces, boundaries, &NCELLS, &NDEGEN_CELLS, &NFACES);
		if (err != 0)
		{
			fprintf(stderr, "build_matrix failed with error code %d\n", err);
			return 1;
		}

		// Solve Linear System A*phi = b for phi
		lapack_int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, A, lda, ipiv, b, ldb);
		if (info != 0)
		{
			fprintf(stderr, "LAPACKE_dgesv failed with error code %d\n", info);
			return 1;
		}
		
		// Lapack dgesv overwrites the right-hand side vector b with the solution, so we can copy it back to phi for the next iteration
		for (int j=0; j < NSOLCELLS; j++)
		{
			phi[j+NDEGEN_CELLS] = b[j];
		}

	}
	//---------------------------------------------------


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
	*			af(owner) += -gamma(face)*Ef(face)/dCF
	*			ac(owner) = gamma(face)*Ef(face)/dCF
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
	free(A);
	free(b);
	free(ipiv);

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

