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

	// Calculate grid dependent coefficient for pressure correction equation stabilization
	double epsilon;
	err = calc_epsilon(cells, &NCELLS, &NDEGEN_CELLS, &epsilon);
	if (err != 0)
	{
		fprintf(stderr, "calc_epsilon failed with error code %d\n", err);
		return 1;
	}

	//double alpha = ALPHA; // velocity urf
	EPSILON_G = epsilon;
	double alpha = ALPHA;
	GAMMA[PCORR] = epsilon + alpha;
	
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

	#if SIMPLE
		double* p = malloc((NCELLS) * sizeof(double)); // pressure array
		if (p == NULL)
		{
			// Print error message to stderr stream and exit
			fprintf(stderr, "Error: Memory allocation failed for pressure array.\n");
			return 1; // Exit with error code
		}

		double* grad_p = malloc((3 * NCELLS) * sizeof(double)); // pressure gradient array
		if (grad_p == NULL)
		{
			// Print error message to stderr stream and exit
			fprintf(stderr, "Error: Memory allocation failed for pressure gradient array.\n");
			return 1; // Exit with error code
		}
		//memset(p, 0, (NCELLS) * sizeof(double)); // Initialize pressure array to zero
		for (int i = 0; i < NCELLS; i++)
		{
			p[i] = 2.0;
		}

		memset(grad_p, 0, (3 * NCELLS) * sizeof(double)); // Initialize pressure gradient array to zero

		
	#endif

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
	boundary boundaries[NBOUNDARIES]; // boundaries

	/* ------------------------ Apply Boundary Conditions ----------------------- */
	for (int i = 0; i < NENTITIES - 1; i++)
	{
		if (cellEntities[i].id != 9) //make entity 9 always internal domain
		{
			// build_boundary_entity(&boundaries[i], i, p1_boundaries[i], p1_boundary_data[i], nodes, faces, cells, &NFACES);
			// Each booundary entity is passed data for every equation
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
#if SIMPLE
			// For SIMPLE, need to pass pressure and pressure gradient to applyBoundary for use in momentum equation boundary conditions
			err = applyBoundary(&boundaries[j], cells, faces, phi[i], grad[i], p, grad_p, GAMMA[i], i, &NCELLS);
#else
			err = applyBoundary(&boundaries[j], cells, faces, phi[i], grad[i], NULL, NULL, GAMMA[i], i, &NCELLS);
#endif // SIMPLE
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
				#if SIMPLE
					// For SIMPLE, need to pass pressure and pressure gradient to applyBoundary for use in
					// momentum equation boundary conditions
					 err = applyBoundary(&boundaries[k], cells, faces, phi[neqn], grad[neqn], p, grad_p, GAMMA[neqn], neqn, &NCELLS);
				#else
					err = applyBoundary(&boundaries[k], cells, faces, phi[neqn], grad[neqn], NULL, NULL, GAMMA[neqn], neqn, &NCELLS);
				#endif // SIMPLE
				if (err != 0)
				{
					fprintf(stderr, "initBoundary failed with error code %d\n", err);
					return 1;
				}
			}
#if SIMPLE
			if (neqn == XMOM)
			{
				// Compute on x  momentum eqn since it wont change for the y momentum
				// For SIMPLE, need to compute pressure gradient before building momentum equation matrix and source term

				// update pressure on boundaries before building momentum equation matrix and source term because pressure gradient is part of the source term for the momentum equation in SIMPLE
				for (int k = 0; k < NBOUNDARIES; k++)
				{
					err = applyBoundary(&boundaries[k],
										cells,
										faces,
										phi[PCORR],
										grad[PCORR],
										p,
										grad_p,
										GAMMA[PCORR],
										PCORR,
										&NCELLS);
					if (err != 0)
					{
						fprintf(stderr, "pressure boundary update failed with error code %d\n", err);
						return 1;
					}
				}
				
				err = compute_lsq_gradient(nodes, cells, faces, &NCELLS, &NDEGEN_CELLS, &NFACES, p, grad_p);
				if (err != 0)
				{
					fprintf(stderr, "compute_lsq_gradient failed with error code %d for pressure gradient\n", err);
					return 1;
				}
			}
#endif //SIMPLE

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

			if (neqn == XMOM || neqn == YMOM)
			{		

				// Pass -GAMMA becasue build diffusion assumes form of the heat equation
				err = build_diffusion(&A[neqn], &b[neqn], phi[neqn], grad[neqn], -GAMMA[neqn], nodes, cells, faces, boundaries, neqn, &NCELLS, &NDEGEN_CELLS, &NFACES);
				if (err != 0)
				{
					fprintf(stderr, "build_diffusion failed with error code %d\n", err);
					return 1;
				}

				err = build_source(&b[neqn], phi[neqn], grad[neqn], neqn, nodes, cells, faces, boundaries, &NCELLS, &NDEGEN_CELLS, &NFACES);
				if (err != 0)
				{
					fprintf(stderr, "build_source failed with error code %d\n", err);
					return 1;
				}

				err = build_gradient(&b[neqn], grad_p, neqn, cells, &NCELLS, &NDEGEN_CELLS);
				if (err != 0)
				{
					fprintf(stderr, "build_gradient failed with error code %d\n", err);
					return 1;
				}
			}
			else if (neqn == PCORR)
			{
				// Recompute Velocity gradient for pressure correction source term since it changes after momentum equation update in SIMPLE
				err = compute_lsq_gradient(nodes, cells, faces,
                               &NCELLS, &NDEGEN_CELLS, &NFACES,
                               phi[XMOM], grad[XMOM]);

    			err = compute_lsq_gradient(nodes, cells, faces,
                               &NCELLS, &NDEGEN_CELLS, &NFACES,
                               phi[YMOM], grad[YMOM]);

				// lhs diffusion of phi_k. Use negative diffusion coefficient because build_diffusion assumes form of the heat equation with -gamma*laplacian(phi) but 
				// pressure correction eqn is gamma*lap(phi)
				err = build_diffusion(&A[neqn], &b[neqn], phi[neqn], grad[neqn], -GAMMA[neqn], nodes, cells, faces, boundaries, neqn, &NCELLS, &NDEGEN_CELLS, &NFACES);
				if (err != 0)
				{
					fprintf(stderr, "build_diffusion failed with error code %d\n", err);
					return 1;
				}

				// rhs div u*
				int eq_ids[ND] = {XMOM, YMOM}; // For divergence, need to pass the velocity equation IDs to build_div to get the correct velocity components
				err = build_div(&b[neqn], phi, grad, eq_ids, nodes, cells, faces, &NCELLS, &NDEGEN_CELLS, &NFACES);
				//err = build_div(&b[neqn], phi, grad, nodes, cells, faces, boundaries, &NCELLS, &NDEGEN_CELLS, &NFACES);
				if (err != 0)
				{
					fprintf(stderr, "build_divergence failed with error code %d\n", err);
					return 1;
				}

				// rhs epsilon * laplacian of p (stabilization term)
				err = build_lap(&b[neqn], p, grad_p, -epsilon, nodes, cells, faces, boundaries, &NCELLS, &NDEGEN_CELLS, &NFACES);
				if (err != 0)
				{
					fprintf(stderr, "build_lap failed with error code %d\n", err);
					return 1;
				}

				// g source term
				err = build_source(&b[neqn], phi[neqn], grad[neqn], neqn, nodes, cells, faces, boundaries, &NCELLS, &NDEGEN_CELLS, &NFACES);
				if (err != 0)
				{
					fprintf(stderr, "build_source failed with error code %d\n", err);
					return 1;
				}
			}
			else
			{
				fprintf(stderr, "Error: Unknown equation index %d\n", neqn);
				return 1;
			}
			/* err = build_diffusion(&A[neqn], &b[neqn], phi[neqn], grad[neqn], GAMMA[neqn], nodes, cells, faces, boundaries, neqn, &NCELLS, &NDEGEN_CELLS, &NFACES);
			if (err != 0)
			{
				fprintf(stderr, "build_diffusion failed with error code %d\n", err);
				return 1;
			}

			err = build_source(&A[neqn], &b[neqn], phi[neqn], grad[neqn],neqn, nodes, cells, faces, boundaries, &NCELLS, &NDEGEN_CELLS, &NFACES);
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
			} */

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
			PetscCall(MatAssemblyBegin(A[neqn], MAT_FINAL_ASSEMBLY));
			PetscCall(MatAssemblyEnd(A[neqn], MAT_FINAL_ASSEMBLY));
			PetscCall(VecAssemblyBegin(b[neqn]));
			PetscCall(VecAssemblyEnd(b[neqn]));

			if (neqn == PCORR)
			{
				PetscInt ref_row = 0; // Index of the row to set as reference for pressure correction (can be any cell index, but typically a boundary cell is chosen)
				PetscScalar ref_value = 0.0; // Reference value for pressure correction

				// Set the diagonal coefficient of the reference row to 1 to specify its value
				// and make pressure unique
				PetscCall(MatZeroRows(A[neqn], 1, &ref_row, 1.0, NULL, NULL));

				PetscCall(VecSetValue(b[neqn], ref_row, ref_value, INSERT_VALUES));
				PetscCall(VecAssemblyBegin(b[neqn]));
				PetscCall(VecAssemblyEnd(b[neqn]));

				// Re-finalize matrix assembly
				PetscCall(MatAssemblyBegin(A[neqn], MAT_FINAL_ASSEMBLY));
				PetscCall(MatAssemblyEnd(A[neqn], MAT_FINAL_ASSEMBLY));
			}

			

			// Set nullspace
			/* MatNullSpace pcorr_nullspace = NULL;
			if (neqn == PCORR)
			{
				PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD,
					PETSC_TRUE,
					0,
					NULL,
					&pcorr_nullspace));
					
					PetscCall(MatSetNullSpace(A[neqn], pcorr_nullspace));

					PetscCall(MatNullSpaceRemove(pcorr_nullspace, b[neqn]));
			} */

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

			/* if (pcorr_nullspace)
			{
				PetscCall(MatNullSpaceDestroy(&pcorr_nullspace));
			} */

			// copy solution back to phi
			const PetscScalar *xarray;
			PetscCall(VecGetArrayRead(xp[neqn], &xarray)); // get pointer to values inside xp vector
			// VecGetArray lets you change the values inside while VecGetArrayRead lets you only read them I think
			for (int k = 0; k < NSOLCELLS; k++)
			{
				phi[neqn][k + NDEGEN_CELLS] = xarray[k];
			}

			PetscCall(VecRestoreArrayRead(xp[neqn], &xarray)); // close access to the values
			PetscCall(KSPDestroy(&ksp[neqn]));

		/* -------------------------------------------------------------------------- */

#endif //(!EXPLICIT && TRANSIENT) || !TRANSIENT

		// Update time if transient
#if TRANSIENT
			time += dt;
#endif // TRANSIENT

		}
		
		// Update solution fields with correction if SIMPLE
#if SIMPLE

		/*
		* Refresh phi_k boundary values after solving PCORR.
		* For homogeneous Neumann, this should set phi_k_boundary = phi_k_owner.
		*/
		for (int k = 0; k < NBOUNDARIES; k++)
		{
			err = applyBoundary(&boundaries[k],
								cells,
								faces,
								phi[PCORR],
								grad[PCORR],
								p,
								grad_p,
								GAMMA[PCORR],
								PCORR,
								&NCELLS);
			if (err != 0)
			{
				fprintf(stderr, "applyBoundary failed for phi_k with error code %d\n", err);
				return 1;
			}
		}

		err = compute_lsq_gradient(nodes, cells, faces,
								&NCELLS, &NDEGEN_CELLS, &NFACES,
								phi[PCORR], grad[PCORR]);
		if (err != 0)
		{
			fprintf(stderr, "compute_lsq_gradient failed for phi_k with error code %d\n", err);
			return 1;
		}

		for (int cell_i = NDEGEN_CELLS; cell_i < NCELLS; cell_i++)
		{
			p[cell_i] += phi[PCORR][cell_i];

			phi[XMOM][cell_i] += -alpha * grad[PCORR][IDX(0, cell_i, 3)];
			phi[YMOM][cell_i] += -alpha * grad[PCORR][IDX(1, cell_i, 3)];
		}

		// Reset pressure correction
		for (int cell_i = 0; cell_i < NCELLS; cell_i++)
		{
			phi[PCORR][cell_i] = 0.0;

			grad[PCORR][IDX(0, cell_i, 3)] = 0.0;
			grad[PCORR][IDX(1, cell_i, 3)] = 0.0;
			grad[PCORR][IDX(2, cell_i, 3)] = 0.0;
		}
#endif //SIMPLE
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
		// update boundaries and gradients before calculating residuals
		for (int k = 0; k < NBOUNDARIES; k++)
		{
			applyBoundary(&boundaries[k], cells, faces,
						phi[XMOM], grad[XMOM],
						p, grad_p,
						GAMMA[XMOM], XMOM, &NCELLS);

			applyBoundary(&boundaries[k], cells, faces,
						phi[YMOM], grad[YMOM],
						p, grad_p,
						GAMMA[YMOM], YMOM, &NCELLS);

			applyBoundary(&boundaries[k], cells, faces,
						phi[PCORR], grad[PCORR],
						p, grad_p,
						GAMMA[PCORR], PCORR, &NCELLS);
		}

		compute_lsq_gradient(nodes, cells, faces,
							&NCELLS, &NDEGEN_CELLS, &NFACES,
							phi[XMOM], grad[XMOM]);

		compute_lsq_gradient(nodes, cells, faces,
							&NCELLS, &NDEGEN_CELLS, &NFACES,
							phi[YMOM], grad[YMOM]);

		compute_lsq_gradient(nodes, cells, faces,
							&NCELLS, &NDEGEN_CELLS, &NFACES,
							p, grad_p);

		double residual[NEQNS];
		bool stop_calc = false;
		bool continue_calc = true;

		/* for (int neqn = 0; neqn < NEQNS; neqn++)
		{
			err = calc_l2_norm(&A[neqn], &b[neqn], phi[neqn], NCELLS, NDEGEN_CELLS, &residual[neqn]);
			//err = calc_Residual(&A[neqn], &b[neqn], phi[neqn], cells, faces, &NCELLS, &NDEGEN_CELLS, &NFACES, &residual[neqn]);
			if (i % RPRT_INTERVAL == 0)
			{
				printf("ITER = %d \n", i + 1);
				printf("Residual = %g \n", residual[neqn]);
			}
		} */
		// Calculate residuals for each eqn
		int eq_ids[ND] = {XMOM, YMOM}; // For divergence, need to pass the velocity equation IDs to build_div to get the correct velocity components

		err = mom_l2_residual(phi[XMOM], grad[XMOM], XMOM, p, grad_p, GAMMA[XMOM], cells, faces, NCELLS, NDEGEN_CELLS, NFACES, &residual[XMOM]);
		err = mom_l2_residual(phi[YMOM], grad[YMOM], YMOM, p, grad_p, GAMMA[YMOM], cells, faces, NCELLS, NDEGEN_CELLS, NFACES, &residual[YMOM]);
		err = continuity_l2_residual(phi, grad, p, grad_p, eq_ids, cells, faces, &NCELLS, &NDEGEN_CELLS, &NFACES, epsilon, &residual[PCORR]);

		for (int neqn = 0; neqn < NEQNS; neqn++)
		{
			if (i % RPRT_INTERVAL == 0)
			{
				printf("ITER = %d \n", i + 1);
				printf("Residual = %g \n", residual[neqn]);
			}
		}

		// if biggest residual is less than stopping condition, stop calculation and print residuals
		if (arr_max(residual, NEQNS) < STOP_COND)
		{
			printf("STOP_COND HIT after %d iterations\n", i + 1);
			for (int neqn = 0; neqn < NEQNS; neqn++)
			{
				printf("Residual of eqn [%d] = %g \n", neqn, residual[neqn]);
			}
			break;
		}
		

	}

	//---------------------------------------------------

	// ------ Write output file --------"output_file.vtk"
#if SIMPLE
	// If using SIMPLE, include pressure and pressure gradient in output for post-processing and visualization purposes
	err = write_vtk_output(out_fname, nodes, cells, &NPOINTS, &NCELLS,
								   &CELL_LIST_SIZE, phi, grad, p, grad_p);
#else
	err = write_vtk_output(out_fname, nodes, cells, &NPOINTS, &NCELLS,
						   &CELL_LIST_SIZE, phi, grad);
#endif //SIMPLE
	if (err != 0)
	{
		fprintf(stderr, "write_vtk_output failed with error code %d\n", err);
		return 1;
	}
#endif // TRANSIENT

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
		//PetscCall(KSPDestroy(&ksp[neqn]));
		PetscCall(VecDestroy(&xp[neqn]));
		PetscCall(VecDestroy(&b[neqn]));
		PetscCall(MatDestroy(&A[neqn]));
	}
	
#if SIMPLE
	free(p);
	free(grad_p);
#endif

	PetscCall(PetscFinalize());
	printf("To C or not to C: that is the question. \n");
	return 0;
}


/*---------------------------------------------------------------------------
* Write data function and grid
* Pass p = NULL and grad_p = NULL if not using SIMPLE
----------------------------------------------------------------------------*/
int write_vtk_output(const char *out_filename, const node *nodes, const cell *cells,
					 const int *NPOINTS, const int *NCELLS, const int *CELL_LIST_SIZE,
					 double **phi, double **grad, const double *p, const double *grad_p)
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
				nodes[i].x, nodes[i].y, nodes[i].z);
	}

	// Write cells
	fprintf(fp, "CELLS %d %d\n", *NCELLS, *CELL_LIST_SIZE);
	for (int i = 0; i < *NCELLS; i++)
	{
		int num_nodes = cells[i].num_nodes;
		fprintf(fp, "%d ", num_nodes);

		for (int point = 0; point < num_nodes; point++)
		{
			if (point == num_nodes - 1)
			{
				fprintf(fp, "%d\n", cells[i].node_ids[point]);
			}
			else
			{
				fprintf(fp, "%d ", cells[i].node_ids[point]);
			}
		}
	}

	// Write cell types
	fprintf(fp, "CELL_TYPES %d\n", *NCELLS);
	for (int i = 0; i < *NCELLS; i++)
	{
		fprintf(fp, "%d\n", cells[i].type);
	}

	// Write cell data for ALL cells, including degenerate cells
	fprintf(fp, "CELL_DATA %d\n", *NCELLS);

	// Cell IDs
	fprintf(fp, "SCALARS cellId int 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (int i = 0; i < *NCELLS; i++)
	{
		fprintf(fp, "%d\n", cells[i].id);
	}

	// Cell Entity Ids
	fprintf(fp, "SCALARS CellEntityIds int 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (int i = 0; i < *NCELLS; i++)
	{
		fprintf(fp, "%d\n", cells[i].entity_id);
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
		fprintf(fp, "%.15f\n", cells[i].volume);
	}

	// Pressure (if SIMPLE)
	if (p != NULL)
	{
		fprintf(fp, "SCALARS pressure double 1\n");
		fprintf(fp, "LOOKUP_TABLE default\n");
		for (int i = 0; i < *NCELLS; i++)
		{
			fprintf(fp, "%.15f\n", p[i]);
		}
	}

	// Pressure gradient
	if (grad_p != NULL)
	{
		// Pressure gradient vector
		fprintf(fp, "VECTORS grad_pressure double\n");
		for (int i = 0; i < *NCELLS; i++)
		{
			fprintf(fp, "%.15f %.15f %.15f\n",
					grad_p[i * 3],
					grad_p[i * 3 + 1],
					grad_p[i * 3 + 2]);
		}
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
		fprintf(fp, "%d\n", (cells[i].volume <= 0.0) ? 1 : 0);
	}

	// Point data
	fprintf(fp, "POINT_DATA %d\n", *NPOINTS);
	fprintf(fp, "SCALARS nodeId int 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (int i = 0; i < *NPOINTS; i++)
	{
		fprintf(fp, "%d\n", nodes[i].id);
	}

	fclose(fp);
	return 0;
}
