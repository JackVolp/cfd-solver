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
#include "setup.h"
#include "math_helpers.h"
#include <petscksp.h>

int main(int argc, char **argv)
{
	/* Setup Petsc */
	PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
	
	// Load grid
	node *nodes;
	cell *cells;
	face *faces;
	cellEntity *cellEntities;

	int NPOINTS = 0, NCELLS = 0, CELL_LIST_SIZE = 0, MAX_FACES = 0, NFACES = 0, NDEGEN_CELLS = 0, NSOLCELLS = 0, NENTITIES = 0;

	/*----------Read grid from .vtk grid file----------*/
	// Load grid from file and store in nodes and cells arrays, also calculate MAX_FACES for memory allocation of faces array
	int err = read_grid(filename, &nodes, &cells, &cellEntities, &NPOINTS, &NCELLS, &CELL_LIST_SIZE, &MAX_FACES, &NDEGEN_CELLS, &NENTITIES);
	NSOLCELLS = NCELLS - NDEGEN_CELLS; // Number of cells that have volume and are included in the solution

	if (err != 0)
	{
		fprintf(stderr, "read_grid failed with error code %d\n", err);
		return 1;
	}

	// Calculate Cell Centroid, Volume, Face information, and other geometric properties
	err = build_faces_and_cells(nodes, cells, &NCELLS, &MAX_FACES, &NFACES, &faces);

	/*----------Allocate Arrays----------*/
	// Allocate conservative scalars
	double* phi[NEQNS];
	double* phi_old[NEQNS];
	double* grad[NEQNS];

	Mat A[NEQNS];		// petsc matrix
	Vec b[NEQNS], xp[NEQNS]; // petsc vectors
	KSP ksp[NEQNS];	// solver object
	PC pc[NEQNS];		// preconditioner object

	for (int i = 0; i < NEQNS; i++)
	{
		phi[i] = malloc((NCELLS) * sizeof(double));
		if (phi[i] == NULL)
		{
			// Print error message to stderr stream and exit
			fprintf(stderr, "Error: Memory allocation failed for phi[%d] array.\n",i);
			return 1; // Exit with error code
		}

		phi_old[i] = malloc((NCELLS) * sizeof(double));
		if (!phi_old[i])
		{
			// Print error message to stderr stream and exit
			fprintf(stderr, "Error: Memory allocation failed for phi_old[%d] array.\n",i);
			return 1; // Exit with error code
		}

		// Allocate gradient array (3 components for x,y,z)
		grad[i] = malloc((3 * NCELLS) * sizeof(double));
		if (grad[i] == NULL)
		{
			// Print error message to stderr stream and exit
			fprintf(stderr, "Error: Memory allocation failed for grad[%d] array.\n",i);
			return 1; // Exit with error code
		}

		// Create PETSc matrix and vectors, and KSP solver context here, and solve the linear system using PETSc. This will require converting the matrix A and vector b into the appropriate PETSc formats (e.g., sparse format for A). You can use the PETSc functions MatSetValues and VecSetValues to populate the matrix and vector, and then call KSPSolve to solve the system. Remember to destroy the PETSc objects after use to free memory.
		// PETSC_COMM_WORLD tells Petsc to use all process in a given run
		// PETSC Decide tells petsc to decide itself how the vector is stored
		// matrix A
		PetscCall(MatCreate(PETSC_COMM_WORLD, &A[i]));
		PetscCall(MatSetSizes(A[i], PETSC_DECIDE, PETSC_DECIDE, NSOLCELLS, NSOLCELLS));
		PetscCall(MatSetFromOptions(A[i]));
		PetscCall(MatSetUp(A[i]));
		// vector b
		PetscCall(VecCreate(PETSC_COMM_WORLD, &b[i]));
		PetscCall(VecSetSizes(b[i], PETSC_DECIDE, NSOLCELLS)); 
		PetscCall(VecSetFromOptions(b[i]));
		PetscCall(VecDuplicate(b[i], &xp[i])); // Create solution vector xp with same size as bp
	}
	/* double *phi = malloc((NEQNS * NCELLS) * sizeof(double));
	if (phi == NULL)
	{
		// Print error message to stderr stream and exit
		fprintf(stderr, "Error: Memory allocation failed for phi array.\n");
		return 1; // Exit with error code
	} */

	/* double *phi_old = malloc((NEQNS * NCELLS) * sizeof(double));
	if (!phi_old)
	{
		// Print error message to stderr stream and exit
		fprintf(stderr, "Error: Memory allocation failed for phi_old array.\n");
		return 1; // Exit with error code
	} */

	// Allocate gradient array (3 components for x,y,z)
	/* double *grad = malloc((3 * NCELLS) * sizeof(double));
	if (grad == NULL)
	{
		// Print error message to stderr stream and exit
		fprintf(stderr, "Error: Memory allocation failed for grad array.\n");
		return 1; // Exit with error code
	} */

	/*--------Setup matrix arrays and solver parameters--------*/
	/* lapack_int n = NSOLCELLS;	// Number of equations (size of the system)
	lapack_int nrhs = 1;		// Number of right-hand sides (columns of B)
	lapack_int lda = NSOLCELLS; // Leading dimension of A
	lapack_int ldb = NSOLCELLS; // Leading dimension of B

	lapack_int *ipiv = malloc(NSOLCELLS * sizeof(lapack_int)); // Pivot indices for LU factorization
	if (!ipiv)
	{
		// Print error message to stderr stream and exit
		fprintf(stderr, "Error: Memory allocation failed for ipiv array.\n");
		return 1; // Exit with error code
	} */

	/* double *A = malloc((NEQNS * NSOLCELLS * NSOLCELLS) * sizeof(double)); // Coefficient matrix (will be stored in sparse format later)
	if (!A)
	{
		// print error message to stderr stream
		fprintf(stderr, "Error: Memory allocation failed for A array.\n");
		return 1; // Exit with error code
	}

	double *b = malloc(NSOLCELLS * sizeof(double)); // Source term vector
	if (!b)
	{
		// print error message to stderr stream
		fprintf(stderr, "Error: Memory allocation failed for b array.\n");
		return 1; // Exit with error code
	} */

	/*--------Initialize Phi--------*/

	for (int i = 0; i < NEQNS; i++)
	{
		// initialize phi to 0 everywhere
		memset(phi[i], 0, (NCELLS) * sizeof(double));

		// initialize grad to zero
		memset(grad[i], 0, ((int)3 * NCELLS) * sizeof(double));
	}
	

/* -------------------------------------------------------------------------- */
/* Initialize Time if transient */
/* -------------------------------------------------------------------------- */
#if TRANSIENT
	double time = 0.0; // Initialize time
	double dt;
	double next_save_time = 0.0; // Initialize next save time for output

#endif

	/*-------- Create and apply boundary conditions--------*/
	// Initialize boundaries (change to allocate for more complex gemoetry)
	// probably move this to setup somehow also
	boundary boundaries[3]; // boundaries

	/* ------------------------ Apply Boundary Conditions ----------------------- */
	for (int i = 0; i < NENTITIES - 1; i++)
	{
		if (cellEntities[i].id != 9) //make entity 9 always internal domain
		{
			// build_boundary_entity(&boundaries[i], i, p1_boundaries[i], p1_boundary_data[i], nodes, faces, cells, &NFACES);
			build_boundary_entity(&boundaries[i], i, problem_boundary_types[i], problem_boundary_data[i], &cellEntities[i], faces, &NFACES);
		}
	}

	//NBOUNDARIES = NENTITIES - 1; // In this case we are treating each entity as a boundary, except for the last one which is the interior cells. This will need to be modified for more complex geometries where not every entity is a boundary
	/* -------------------------------------------------------------------------- */
	/* Solver Loop */
	/* -------------------------------------------------------------------------- */
	printf("Start Solving \n");

	// Apply boundary conditions (sets phi on boundaries)
	for (int i = 0; i < NEQNS; i++)
	{
		for (int j = 0; j < NBOUNDARIES; j++)
		{
			err = applyBoundary(&boundaries[j], cells, faces, phi[i], grad[i], GAMMA[i], &NCELLS);

			if (err != 0)
			{
				fprintf(stderr, "initBoundary failed with error code %d\n", err);
				return 1;
			}
		}
	}
	

	for (int i = 0; i < MAX_ITER; i++)
	{
		for (int neqn = 0; neqn < NEQNS; neqn++)
		{
			// Save old phi
			phi_old[neqn] = memcpy(phi_old[neqn], phi[neqn], (NCELLS) * sizeof(double));

			// Apply boundary conditions (sets phi on boundaries)
			for (int k = 0; k < NBOUNDARIES; k++)
			{
				err = applyBoundary(&boundaries[k], cells, faces, phi[neqn], grad[neqn], GAMMA[neqn], &NCELLS);

				if (err != 0)
				{
					fprintf(stderr, "initBoundary failed with error code %d\n", err);
					return 1;
				}
			}

			// Compute Gradient at cell centers and at boundary faces
			err = compute_lsq_gradient(nodes, cells, faces, &NCELLS, &NDEGEN_CELLS, &NFACES, phi[neqn], grad[neqn]);
			if (err != 0)
			{
				fprintf(stderr, "compute_lsq_gradient failed with error code %d for equation %d\n", err, neqn);
				return 1;
			}

			// Build Matrix and Source term
			// initialize matrix coefficients and source term vector to zero
			PetscCall(MatZeroEntries(A[neqn]));
			PetscCall(VecZeroEntries(b[neqn]));
		
			/* memset(A, 0, (NEQNS * NSOLCELLS * NSOLCELLS) * sizeof(double));
			memset(b, 0, ((NEQNS * NSOLCELLS) * sizeof(double))); */

			err = build_diffusion(&A[neqn], &b[neqn], phi[neqn], grad[neqn], GAMMA[neqn], nodes, cells, faces, boundaries, &NCELLS, &NDEGEN_CELLS, &NFACES);
			if (err != 0)
			{
				fprintf(stderr, "build_diffusion failed with error code %d\n", err);
				return 1;
			}

			err = build_source(&A[neqn], &b[neqn], phi[neqn], grad[neqn], nodes, cells, faces, boundaries, &NCELLS, &NDEGEN_CELLS, &NFACES);
			if (err != 0)
			{
				fprintf(stderr, "build_source failed with error code %d\n", err);
				return 1;
			}

			err = build_advection(&A[neqn], &b[neqn], phi[neqn], grad[neqn], nodes, cells, faces, boundaries, &NCELLS, &NDEGEN_CELLS, &NFACES);
			if (err != 0)
			{
				fprintf(stderr, "build_advection failed wiht error code %d\n", err);
				return 1;
			}

#if TRANSIENT

			// Intermediate matrix assembly so matrix A can be indexed.
			PetscCall(MatAssemblyBegin(A[neqn], MAT_FLUSH_ASSEMBLY));
			PetscCall(MatAssemblyEnd(A[neqn], MAT_FLUSH_ASSEMBLY));

			// Calculate time step here
			err = calc_time_step(cells, &A[neqn], &NCELLS, &NDEGEN_CELLS, &time, &dt);
			if (err != 0)
			{
				fprintf(stderr, "calc_time_step failed with error code %d\n", err);
				return 1;
			}

#if EXPLICIT
			err = explicit_update(&A[i], &b[i], phi[i], cells, faces, &dt, &NCELLS, &NDEGEN_CELLS);
			if (err != 0)
			{
				fprintf(stderr, "explicit_update failed with error code %d\n", err);
				return 1;
			}
#else
			// build transient contribution to matrix and source term
			err = build_transient(&A[neqn], &b[neqn], phi[neqn], cells, &NCELLS, &NDEGEN_CELLS, dt);
			if (err != 0)
			{
				fprintf(stderr, "build_transient failed with error code %d\n", err);
				return 1;
			}
#endif // EXPLICIT
#endif // TRANSIENT

#if (!EXPLICIT && TRANSIENT) || !TRANSIENT
			/* -------------------------------------------------------------------------- */
			/*                       Solve Linear System                      */
			/* -------------------------------------------------------------------------- */
			/* --------------------------------- LAPACK --------------------------------- */
			// Explicit update of phi using forward Euler time stepping

			// Solve Linear System A*phi = b for phi
			/* lapack_int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, A, lda, ipiv, b, ldb);
			if (info != 0)
			{
				fprintf(stderr, "LAPACKE_dgesv failed with error code %d\n", info);
				return 1;
			}

			// Lapack dgesv overwrites the right-hand side vector b with the solution, so we can copy it back to phi for the next iteration. We can also use b to compute the maximum change
			for (int j = 0; j < NSOLCELLS; j++)
			{
				phi[j + NDEGEN_CELLS] = b[j];
			} */

			/* ---------------------------------- PETSc --------------------------------- */

			/* Mat Ap;		// petsc matrix
			Vec bp, xp; // petsc vectors
			KSP ksp;	// solver object
			PC pc;		// preconditioner object

			// Create PETSc matrix and vectors, and KSP solver context here, and solve the linear system using PETSc. This will require converting the matrix A and vector b into the appropriate PETSc formats (e.g., sparse format for A). You can use the PETSc functions MatSetValues and VecSetValues to populate the matrix and vector, and then call KSPSolve to solve the system. Remember to destroy the PETSc objects after use to free memory.
			// matrix A
			PetscCall(MatCreate(PETSC_COMM_WORLD, &Ap));
			PetscCall(MatSetSizes(Ap, PETSC_DECIDE, PETSC_DECIDE, NSOLCELLS, NSOLCELLS));
			PetscCall(MatSetFromOptions(Ap));
			PetscCall(MatSetUp(Ap));
			// vector b
			PetscCall(VecCreate(PETSC_COMM_WORLD, &bp));
			PetscCall(VecSetSizes(bp, PETSC_DECIDE, NSOLCELLS));
			PetscCall(VecSetFromOptions(bp));
			PetscCall(VecDuplicate(bp, &xp)); // Create solution vector xp with same size as bp
	*/

			// Copy dense matrix into petsc matrix
			/* for (int row = 0; row < NSOLCELLS; row++)
			{
				PetscCall(VecSetValue(b, row, b[row], INSERT_VALUES));

				for (int col = 0; col < NSOLCELLS; col++)
				{
					double value = A[row + col * NSOLCELLS]; // LAPACK_COL_MAJOR layout 

					if (fabs(value) > 0.0)
					{
						PetscCall(MatSetValue(Ap, row, col, value, INSERT_VALUES));
					}
				}
			} */

			PetscCall(MatAssemblyBegin(A[neqn], MAT_FINAL_ASSEMBLY));
			PetscCall(MatAssemblyEnd(A[neqn], MAT_FINAL_ASSEMBLY));
			PetscCall(VecAssemblyBegin(b[neqn]));
			PetscCall(VecAssemblyEnd(b[neqn]));

			// Solve system
			PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp[neqn]));
			PetscCall(KSPSetOperators(ksp[neqn], A[neqn], A[neqn]));

			if (!PCMPIServerActive)
			{
				PetscCall(KSPGetPC(ksp[neqn], &pc[neqn]));
				PetscCall(PCSetType(pc[neqn], PCILU));	  // set ilu as default preconditioner
				PetscCall(KSPSetType(ksp[neqn], KSPGMRES)); // default solve with gmres
			}
			

			// allow for settings override at runtime
			PetscCall(KSPSetFromOptions(ksp[neqn]));

			PetscCall(KSPSolve(ksp[neqn], b[neqn], xp[neqn])); // solve

			// copy solution back to phi
			const PetscScalar *xarray;
			PetscCall(VecGetArrayRead(xp[neqn], &xarray)); // get pointer to values inside xp vector
			// VecGetArray lets you change the values inside while VecGetArrayRead lets you only read them I think
			for (int k = 0; k < NSOLCELLS; k++)
			{
				phi[neqn][k + NDEGEN_CELLS] = xarray[k];
			}

			PetscCall(VecRestoreArrayRead(xp[neqn], &xarray)); // close access to the values

		/* -------------------------------------------------------------------------- */

#endif //(!EXPLICIT && TRANSIENT) || !TRANSIENT

		// Update time if transient
#if TRANSIENT
			time += dt;
#endif // TRANSIENT

		}
		

		// Stopping conditions
#if TRANSIENT

		// Reporting to console
		if (i % RPRT_INTERVAL == 0)
		{
			printf("ITER = %d \n", i + 1);
			printf("Time = %g \n", time);
		}

		// Saving solution based on TIME_INTERVAL
		if (time >= next_save_time - 1e-6)
		{
			// const char* out_fname = "hw2_20x20_out.vtk";
			char base[256];
			char out_fname_time[256];

			// Copy string into base
			snprintf(base, sizeof(base), "%s", out_fname);

			// Remove '.vtk' from filename
			char *p_dot = strrchr(base, '.');

			if (p_dot)
			{
				*p_dot = '\0';
			}

			snprintf(out_fname_time, sizeof(out_fname_time), "%s_%04d.vtk", base, i);

			printf("Saving output at time %g to file: %s\n", time, out_fname_time);

			err = write_vtk_output(out_fname_time, &nodes, &cells, &NPOINTS, &NCELLS,
								   &CELL_LIST_SIZE, phi, grad);
			if (err != 0)
			{
				fprintf(stderr, "write_vtk_output failed with error code %d\n", err);
				return 1;
			}
			next_save_time += SAVE_INTERVAL; // Update next save time
		}
		if (time >= T_FINAL)
			break;
#else
		// Residual/linear system based stopping condition.
		double residual[NEQNS];
		bool stop_calc = false;
		bool continue_calc = true;

		for (int neqn = 0; neqns < NEQNS; neqns++)
		{
			err = calc_Residual(&A, &b, phi, cells, faces, &NCELLS, &NDEGEN_CELLS, &NFACES, &residual);
			if (i % RPRT_INTERVAL == 0)
			{
				printf("ITER = %d \n", i + 1);
				printf("Residual = %g \n", residual);
			}
		}
		
		if (all_greater_than_or_eq(residual, NEQNS, STOP_COND))
		{
			printf("STOP_COND HIT after %d iterations\n", i + 1);
			for (int neqn = 0; neqn < NEQNS; neqn++)
			{
				printf("Residual of eqn [%d] = %g \n", neqn, residual[neqn]);
			}
			break;
		}
		
#endif // TRANSIENT
	}

	//---------------------------------------------------

	// ------ Write output file --------"output_file.vtk"
	err = write_vtk_output(out_fname, &nodes, &cells, &NPOINTS, &NCELLS,
						   &CELL_LIST_SIZE, phi, grad);
	if (err != 0)
	{
		fprintf(stderr, "write_vtk_output failed with error code %d\n", err);
		return 1;
	}

	// Release Allocated Memory for grid
	free_grid(nodes, cells, faces, NCELLS, NFACES);

	// Release conservatiWve scalars memory
	for (int neqn = 0; neqn < NEQNS; neqn++)
	{
		free(phi[neqn]);
		free(grad[neqn]);
		//free(A);
		//free(b);
		//free(ipiv);
		free(phi_old[neqn]);

		// Clean up PETSc objects
		PetscCall(KSPDestroy(&ksp[neqn]));
		PetscCall(VecDestroy(&xp[neqn]));
		PetscCall(VecDestroy(&b[neqn]));
		PetscCall(MatDestroy(&A[neqn]));
	}
	

	PetscCall(PetscFinalize());
	printf("To C or not to C: that is the question. \n");
	return 0;
}


/*---------------------------------------------------------------------------
* Write data function and grid
----------------------------------------------------------------------------*/
int write_vtk_output(const char *out_filename, node **nodes, cell **cells,
					 int *NPOINTS, int *NCELLS, int *CELL_LIST_SIZE, double **phi, double **grad)
{
	// Open file for writing
	FILE *fp = fopen(out_filename, "w");

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

	// Cell Entity Ids
	fprintf(fp, "SCALARS CellEntityIds int 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (int i = 0; i < *NCELLS; i++)
	{
		fprintf(fp, "%d\n", (*cells)[i].entity_id);
	}

	// Phi
	for (int neqn = 0; neqn < NEQNS; neqn++)
	{
		fprintf(fp, "SCALARS phi[%d] double 1\n", neqn);
		fprintf(fp, "LOOKUP_TABLE default\n");
		for (int i = 0; i < *NCELLS; i++)
		{
			//fprintf(fp, "%.15f\n", (*(phi+neqn))[IDX(i, 0, *NCELLS)]);
			fprintf(fp, "%.15f\n", phi[neqn][i]);
		}
	}
	

	// Cell volume
	fprintf(fp, "SCALARS cellVolume double 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (int i = 0; i < *NCELLS; i++)
	{
		fprintf(fp, "%.15f\n", (*cells)[i].volume);
	}

	// Gradient vector
	for (int neqn = 0; neqn < NEQNS; neqn++)
	{
		fprintf(fp, "VECTORS grad_phi[%d] double\n",neqn);
		for (int i = 0; i < *NCELLS; i++)
		{
			fprintf(fp, "%.15f %.15f %.15f\n",
					grad[neqn][IDX(0, i, 3)],
					grad[neqn][IDX(1, i, 3)],
					grad[neqn][IDX(2, i, 3)]);
		}
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
