
#include "solver.h"

int compute_lsq_gradient(node* nodes, cell* cells, face* faces, int* NCELLS,
	int* NDEGEN_CELLS, int* NFACES, double* phi, double* grad)

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
		b[IDX(vC_idx,0, NVOL_CELLS)] += w * dphi * dxk; // x component/row of b
		b[IDX(vC_idx,1, NVOL_CELLS)] += w * dphi * dyk; // y component/row of b

		// Update A11, A12, A22, and b for neighbor cell if the neighbor cell is not degenerate. Neighbor will be degenerate for boundary cells.
		if (!f->boundary_face)
		{
			A11[vF_idx] += w * dxk * dxk;
			A12[vF_idx] += w * dxk * dyk; // same as A21
			A22[vF_idx] += w * dyk * dyk;

			// Update b for neighbor cell
			b[IDX(vF_idx, 0, NVOL_CELLS)] += w * dphi * dxk; // x component/row of b
			b[IDX(vF_idx, 1, NVOL_CELLS)] += w * dphi * dyk; // y component/row of b
		}



	}

	// Loop over all cells and solve for gradient
	for (int i = 0; i < NVOL_CELLS; i++)
	{
		int cell_id = i + *NDEGEN_CELLS; // Adjust index to account for degenerate cells at the beginning of the cells array

		solve_2x2_system(A11[i], A12[i], A12[i], A22[i],
			b[IDX(i, 0, NVOL_CELLS)], b[IDX(i, 1, NVOL_CELLS)],
			&grad[IDX(0, cell_id, 3)], &grad[IDX(1, cell_id, 3)]); // Store gradient in correct location in grad array based on cell id
		grad[IDX(2, cell_id, 3)] = 0.0;
	}

	// Free allocated memory for gradient coefficient matrices
	free(A11);
	free(A12);
	free(A22);
	free(b);
	return 0;
}


// this function also updates the gradient vector at the boundary /degenerate cell indicies with the gradients at the boundary faces. This should prob be in a different function
int build_matrix(double* A, double* b, double* phi, double* grad, node* nodes, cell* cells, face* faces, boundary* boundaries, int* NCELLS, int* NDEGEN_CELLS, int* NFACES)
{
	int NSOLCELLS = (*NCELLS) - (*NDEGEN_CELLS); // Number of cells included in solution (non-degenerate cells)
		
	// Loop over all all faces and add matrix contributions
	for (int i = 0; i < *NFACES; i++)
	{
		face* f = &faces[i]; //current face
		cell* cell_C = &cells[f->owner]; //Face owner cell
		cell* cell_F = &cells[f->neighbor]; //Face neighbor cell

		int C_idx = cell_C->id; // Index for owner cell in A and b arrays
		int F_idx = cell_F->id; // Index for neighbor cell in A and b arrays
		
		int Csol_idx = C_idx - *NDEGEN_CELLS; // Index for owner cell in phi and grad arrays (only volume cells are included in solution)

		int Fsol_idx = F_idx - *NDEGEN_CELLS; // Index for neighbor cell in phi and grad arrays (only volume cells are included in solution)

		// Calculate gDiff
		double rCF[3] = { cell_F->xc - cell_C->xc, cell_F->yc - cell_C->yc, cell_F->zc - cell_C->zc };
		double dCF = 0;
		magnitude(rCF, &dCF);

		double Ef[3] = { f->Ex, f->Ey, f->Ez }; // orthogonal-like Surface area vector of face
		double Ef_mag = 0;
		magnitude(Ef, &Ef_mag);

		double gDiff = Ef_mag / dCF; // "Geometric Diffusion Coefficient"

		double Tf[3] = { f->Tx, f->Ty, f->Tz }; // Tangential contribution vector of face

		if (f->boundary_face)
		{
			// make new variables to match book, same as above though smeel
			double dCb = dCF; // Distance from cell centroid to boundary face centroid, should be positive since face vector points outward from owner cell

			double gDiff_b = Ef_mag / dCb; // Geometric diffusion coefficient for boundary face

			// Face area magnitude
			double Sf_mag = mag(((double[3]){ f->Sx, f->Sy, f->Sz }));
 // Magnitude of surface area vector of face

			switch (boundaries[f->boundary_id].type)
			{
				case Dirichlet: {
					
					A[IDX(Csol_idx, Csol_idx, NSOLCELLS)] += GAMMA * gDiff_b; //aC for owner
					
					//interpolate gradient to face using Eq. (9.33) but with boundary value instead of neighbor cell value
					double grad_face[3] = { 0., 0., 0. }; // Initialize gradient at face
					grad2face(grad_face, &grad[3*C_idx], &grad[3*F_idx], rCF, dCF, phi[C_idx], phi[F_idx]);

					// Update Gradient at face
					grad[IDX(0, F_idx, 3)] = grad_face[0];
					grad[IDX(1, F_idx, 3)] = grad_face[1];
					grad[IDX(2, F_idx, 3)] = grad_face[2];

					double fluxVb = -GAMMA * gDiff_b * phi[F_idx] - GAMMA*dot(grad_face,Tf); //nonlinearized flux contibution

					b[Csol_idx] += -fluxVb; // Source term contribution for owner cell, negative since we are moving it to the right hand side of the equation
					break;
				}
				case Neumann: {
					// Nothing needs to be done to coefficients for neumann as long as the source term has already been initialized prior. I.e, b has Q_c*V_c added to it already.
					double q_b = boundaries[f->boundary_id].data.q_b; // Neumann boundary condition value (flux)									

					double fluxVb = q_b * Sf_mag; // Flux contribution from boundary condition, positive since we are adding it to the source term on the right hand side of the equation

					b[Csol_idx] += fluxVb; // Source term contribution for owner cell
					break;
				}
				case Robin: {
					double h_inf = boundaries[f->boundary_id].data.robin.h_inf;
					double phi_inf = boundaries[f->boundary_id].data.robin.phi_inf;

					double fluxCb = (h_inf * Sf_mag * GAMMA * gDiff_b)
						/ (h_inf * Sf_mag + GAMMA * gDiff_b); // Coefficient for phi at owner cell in Robin boundary condition eq. 8.87

					//interpolate gradient to face using Eq. (9.33) but with boundary value instead of neighbor cell value
					double grad_face[3] = { 0., 0., 0. }; // Initialize gradient at face
					grad2face(grad_face, &grad[3 * C_idx], &grad[3 * F_idx], rCF, dCF, phi[C_idx], phi[F_idx]);

					// Update Gradient at face
					grad[IDX(0, F_idx, 3)] = grad_face[0];
					grad[IDX(1, F_idx, 3)] = grad_face[1];
					grad[IDX(2, F_idx, 3)] = grad_face[2];

					double fluxVb = -fluxCb * phi_inf - (h_inf * Sf_mag * GAMMA * dot(grad_face, Tf)) / (h_inf * Sf_mag + GAMMA * gDiff_b); // eq. 8.87, nonlinearized flux contribution 

					A[IDX(Csol_idx, Csol_idx, NSOLCELLS)] += fluxCb; //aC for owner cell from Robin boundary condition
					b[Csol_idx] += -fluxVb; // Source term contribution for owner cell, negative since we are moving it to the right hand side of the equation
					break;
				}
			}
		}
		else // interoir face 
		{					
			// Assign contribution to coefficients for owner cell (cell_C)
			A[IDX(Csol_idx, Csol_idx, NSOLCELLS)] += GAMMA*gDiff; //aC for owner
			A[IDX(Csol_idx, Fsol_idx, NSOLCELLS)] += -GAMMA*gDiff; //aF for owner

			// Assign contribution to coefficients for neighbor cell (cell_F)
			A[IDX(Fsol_idx, Fsol_idx, NSOLCELLS)] += GAMMA*gDiff; //aC for neighbor
			A[IDX(Fsol_idx, Csol_idx, NSOLCELLS)] += -GAMMA*gDiff; //aF for neighbor

			// ----Source Terms------
			// interpolate gradient to face using Eq. (9.33)

			double grad_face[3] = { 0., 0., 0. }; // Initialize gradient at face
			grad2face(grad_face, &grad[3*C_idx], &grad[3*F_idx], rCF, dCF, phi[C_idx], phi[F_idx]);

			// Contribution to source term for owner cell
			
			b[Csol_idx] += GAMMA * dot(grad_face, Tf); // Source term contribution for owner cell
			b[Fsol_idx] += -GAMMA * dot(grad_face, Tf); // Source term contribution for neighbor cell (negative of owner contribution)
		}
	}

	//initialize b with the energy source term for all cells
	for (int i = 0; i < NSOLCELLS; i++)
	{
		cell* c = &cells[i+(*NDEGEN_CELLS)];

		double Q = Q_C(c->xc, c->yc, c->zc);

		b[i] += Q * c->volume;
	}

	return 0;
}

int grad2face(double* grad_face, double* grad_C, double* grad_F, double* rCF, double dCF, double phi_C, double phi_F)
{
	// Compute average gradient at face
	double grad_face_avg[3] = {
		(grad_C[0] + grad_F[0]) / 2.0,
		(grad_C[1] + grad_F[1]) / 2.0,
		(grad_C[2] + grad_F[2]) / 2.0
	}; // Average gradient at face
	double eCF[3] = { rCF[0] / dCF, rCF[1] / dCF, rCF[2] / dCF }; // Unit vector from cell C to cell F
	// Gradient correction
	double correction = (phi_F - phi_C) / dCF - dot(grad_face_avg, eCF); // Average gradient at face dotted with unit vector from cell C to cell F
	double correction_vec[3] = {
		correction * eCF[0], correction * eCF[1], correction * eCF[2]
	}; // Correction vector to be added to gradient at face
	// Final gradient at face after correction
	grad_face[0] = grad_face_avg[0] + correction_vec[0];
	grad_face[1] = grad_face_avg[1] + correction_vec[1];
	grad_face[2] = grad_face_avg[2] + correction_vec[2];
	return 0;
}

int applyBoundary(boundary* b, cell* cells,
	face* faces, double* phi, double* grad, int* NCELLS)
{
	// Loop over all faces in boundary and apply boundary conditions
	for (int i = 0; i < b->num_faces; i++)
	{
		int face_id = b->face_ids[i];
		face* f = &faces[face_id];
		cell* c_owner = &cells[f->owner];

		// Should always be a boundary face but just to be sure
		int phi_face_idx = f->neighbor; // If boundary face, use owner cell for phi index, otherwise use neighbor cell (should not trigger for internal faces)

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
			double gradPhi_face[3] = { grad[IDX(0,phi_face_idx,3)], grad[IDX(1,phi_face_idx,3)], grad[IDX(2,phi_face_idx,3)] };

			double grad_dot_Tf = gradPhi_face[0] * f->Tx + gradPhi_face[1] * f->Ty + gradPhi_face[2] * f->Tz; // Gradient at face dot tangential contribution vector of face

			phi[phi_face_idx] = (h_inf * Sf_mag * phi_inf
				+ GAMMA * gDiff * phi[phi_owner_idx]
				- GAMMA * grad_dot_Tf) / (h_inf * Sf_mag + GAMMA * gDiff);
			break;
		}
		default:
			fprintf(stderr, "Error: Unknown boundary condition type for boundary ID %d\n", b->id);
			return 1;
		}
	}


	return 0;
}
