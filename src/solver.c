
#include "solver.h"
#include "setup.h"
#include <math.h>
#include <string.h>

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
		double dxk = F->xc - C->xc;
		double dyk = F->yc - C->yc;
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


/* -------------------------------------------------------------------------- */
/*              Diffusion term matrix contribution \nabla^2 \phi              */
/* -------------------------------------------------------------------------- */
// this function also updates the gradient vector at the boundary /degenerate cell indicies with the gradients at the boundary faces. This should prob be in a different function
int build_diffusion(Mat* A, Vec* b, double* phi, double* grad, node* nodes, cell* cells, face* faces, boundary* boundaries, int* NCELLS, int* NDEGEN_CELLS, int* NFACES)
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

			boundary bound = boundaries[f->boundary_id];

			switch (bound.type)
			{
				case Dirichlet: {
					
					//A[IDX(Csol_idx, Csol_idx, NSOLCELLS)] += GAMMA * gDiff_b; //aC for owner
					PetscScalar val = GAMMA * gDiff_b; //aC for owner
					PetscCall(MatSetValues(*A, 1, &Csol_idx, 1, &Csol_idx, &val, ADD_VALUES));
					
					//interpolate gradient to face using Eq. (9.33) but with boundary value instead of neighbor cell value
					double grad_face[3] = { 0., 0., 0. }; // Initialize gradient at face
					grad2face(grad_face, &grad[3*C_idx], &grad[3*F_idx], rCF, dCF, phi[C_idx], phi[F_idx],cell_C,cell_F,f);

					// Update Gradient at face
					grad[IDX(0, F_idx, 3)] = grad_face[0];
					grad[IDX(1, F_idx, 3)] = grad_face[1];
					grad[IDX(2, F_idx, 3)] = grad_face[2];

					double fluxVb = -GAMMA * gDiff_b * phi[F_idx] - GAMMA*dot(grad_face,Tf); //nonlinearized flux contibution
					
					val = -fluxVb; //petsc scalar to add to matrix source
					PetscCall(VecSetValues(*b, 1, &Csol_idx, &val, ADD_VALUES));
					//b[Csol_idx] += -fluxVb; // Source term contribution for owner cell, negative since we are moving it to the right hand side of the equation
					break;
				}
				case Neumann: {
					// Nothing needs to be done to coefficients for neumann as long as the source term has already been initialized prior. I.e, b has Q_c*V_c added to it already.

					double q_b = bound.data.q_b(&bound, f, 0.0);
					//double q_b = boundaries[f->boundary_id].data.q_b; // Neumann boundary condition value (flux)									

					double fluxVb = q_b * Sf_mag; // Flux contribution from boundary condition, positive since we are adding it to the source term on the right hand side of the equation
					
					PetscScalar val = fluxVb; //petsc scalar to add to matrix source
					PetscCall(VecSetValues(*b, 1, &Csol_idx, &val, ADD_VALUES));
					//b[Csol_idx] += fluxVb; // Source term contribution for owner cell
					break;
				}
				case Robin: {
					double h_inf = bound.data.robin.h_inf(&bound, f, 0.0); 
					double phi_inf = bound.data.robin.phi_inf(&bound, f, 0.0);

					double fluxCb = (h_inf * Sf_mag * GAMMA * gDiff_b)
						/ (h_inf * Sf_mag + GAMMA * gDiff_b); // Coefficient for phi at owner cell in Robin boundary condition eq. 8.87

					//interpolate gradient to face using Eq. (9.33) but with boundary value instead of neighbor cell value
					double grad_face[3] = { 0., 0., 0. }; // Initialize gradient at face
					grad2face(grad_face, &grad[3 * C_idx], &grad[3 * F_idx], rCF, dCF, phi[C_idx], phi[F_idx],cell_C,cell_F,f);

					// Update Gradient at face
					grad[IDX(0, F_idx, 3)] = grad_face[0];
					grad[IDX(1, F_idx, 3)] = grad_face[1];
					grad[IDX(2, F_idx, 3)] = grad_face[2];

					double fluxVb = -fluxCb * phi_inf - (h_inf * Sf_mag * GAMMA * dot(grad_face, Tf)) / (h_inf * Sf_mag + GAMMA * gDiff_b); // eq. 8.87, nonlinearized flux contribution 
					
					PetscScalar val = fluxCb; //petsc scalar to add to aC
					PetscCall(MatSetValues(*A, 1, &Csol_idx, 1, &Csol_idx, &val, ADD_VALUES));

					val = -fluxVb; //petsc scalar to add to RHS
					PetscCall(VecSetValues(*b, 1, &Csol_idx, &val, ADD_VALUES));
					
					//A[IDX(Csol_idx, Csol_idx, NSOLCELLS)] += fluxCb; //aC for owner cell from Robin boundary condition
					//b[Csol_idx] += -fluxVb; // Source term contribution for owner cell, negative since we are moving it to the right hand side of the equation
					break;
				}
			}
		}
		else // interoir face 
		{					
			// Assign contribution to coefficients for owner cell (cell_C)
			PetscScalar val_aCC = GAMMA*gDiff;
			PetscScalar val_aCF = -GAMMA*gDiff;

			PetscCall(MatSetValues(*A, 1, &Csol_idx, 1, &Csol_idx, &val_aCC, ADD_VALUES));
			PetscCall(MatSetValues(*A, 1, &Csol_idx, 1, &Fsol_idx, &val_aCF, ADD_VALUES));

			//A[IDX(Csol_idx, Csol_idx, NSOLCELLS)] += GAMMA*gDiff; //aC for owner
			//A[IDX(Csol_idx, Fsol_idx, NSOLCELLS)] += -GAMMA*gDiff; //aF for owner

			// Assign contribution to coefficients for neighbor cell (cell_F)
			PetscScalar val_aFF = GAMMA*gDiff;
			PetscScalar val_aFC = -GAMMA*gDiff;

			PetscCall(MatSetValues(*A, 1, &Fsol_idx, 1, &Fsol_idx, &val_aFF, ADD_VALUES));
			PetscCall(MatSetValues(*A, 1, &Fsol_idx, 1, &Csol_idx, &val_aFC, ADD_VALUES));
			//A[IDX(Fsol_idx, Fsol_idx, NSOLCELLS)] += GAMMA*gDiff; //aC for neighbor
			//A[IDX(Fsol_idx, Csol_idx, NSOLCELLS)] += -GAMMA*gDiff; //aF for neighbor

			// ----Source Terms------
			// interpolate gradient to face using Eq. (9.33)

			double grad_face[3] = { 0., 0., 0. }; // Initialize gradient at face
			grad2face(grad_face, &grad[3*C_idx], &grad[3*F_idx], rCF, dCF, phi[C_idx], phi[F_idx],cell_C,cell_F,f);

			// Contribution to source term for owner cell
			PetscScalar val_bC = GAMMA * dot(grad_face, Tf);
			PetscScalar val_bF = -GAMMA * dot(grad_face, Tf);

			PetscCall(VecSetValues(*b, 1, &Csol_idx, &val_bC, ADD_VALUES));
			PetscCall(VecSetValues(*b, 1, &Fsol_idx, &val_bF, ADD_VALUES));

			//b[Csol_idx] += GAMMA * dot(grad_face, Tf); // Source term contribution for owner cell
			//b[Fsol_idx] += -GAMMA * dot(grad_face, Tf); // Source term contribution for neighbor cell (negative of owner contribution)
		}
	}
	return 0;
}

/* -------------------------------------------------------------------------- */
/*                   Soure term matrix contribution routine                   */
/* -------------------------------------------------------------------------- */
int build_source(Mat* A, Vec* b, double* phi, double* grad, node* nodes, cell* cells, face* faces, boundary* boundaries, int* NCELLS, int* NDEGEN_CELLS, int* NFACES)
{
	int NSOLCELLS = (*NCELLS) - (*NDEGEN_CELLS); // Number of cells included in solution (non-degenerate cells)

	for (int i = 0; i < NSOLCELLS; i++)
	{
		cell* c = &cells[i+(*NDEGEN_CELLS)];

		double Q = Q_C(c->xc, c->yc, c->zc);

		PetscScalar val_bC = Q * c->volume;
		PetscCall(VecSetValues(*b, 1, &i, &val_bC, ADD_VALUES));
		//b[i] += Q * c->volume;
	}

	return 0;
}

/* -------------------------------------------------------------------------- */
/*          Routine to add advection term contribution div(rho*u*phi)         */
/* -------------------------------------------------------------------------- */
int build_advection(Mat* A, Vec* b, double* phi, double* grad, node* nodes, cell* cells, face* faces, boundary* boundaries, int* NCELLS, int* NDEGEN_CELLS, int* NFACES)
{
	int NSOLCELLS = (*NCELLS) - (*NDEGEN_CELLS); // Number of cells included in solution (non-degenerate cells)

	//Loop over all faces
	for (int i = 0; i < *NFACES; i++)
	{
		face* f = &faces[i];
		cell* cell_C = &cells[f->owner];
		cell* cell_F = &cells[f->neighbor];

		int C_idx = cell_C->id; // Index for owner cell in A and b arrays
		int F_idx = cell_F->id; // Index for neighbor cell in A and b arrays

		int Csol_idx = C_idx - *NDEGEN_CELLS; // Index for owner cell in phi and grad arrays (only volume cells are included in solution)

		int Fsol_idx = F_idx - *NDEGEN_CELLS; // Index for neighbor cell in phi and grad arrays (only volume cells are included in solution)

		/* -------------------------------------------------------------------------- */
		/* Calculate Face mass flow rate */
		/* -------------------------------------------------------------------------- */
		// Positive is mdot_f is in direction of Sf
		double v[3] = { XVEL, YVEL, 0.0 }; //velocity vector
		double Sf[3] = { f->Sx, f->Sy, f->Sz };

		double mdot_f = RHO * dot(v, Sf); //Mass flow through face
		double Sf_mag = mag(Sf);

		if (f->boundary_face)
		{
			// For outflow case (mdotf will be positive), normal contribution
			PetscScalar val_aCC = fmax(mdot_f, 0.0);
			PetscCall(MatSetValues(*A, 1, &Csol_idx, 1, &Csol_idx, &val_aCC, ADD_VALUES));

			//A[IDX(Csol_idx, Csol_idx, NSOLCELLS)] += fmax(mdot_f, 0.0);

			//inflow case (mdotf is negative), boundary flux added to RHS
			PetscScalar val_bC = fmax(-mdot_f, 0) * phi[F_idx];
			PetscCall(VecSetValues(*b, 1, &Csol_idx, &val_bC, ADD_VALUES));

			//b[Csol_idx] += fmax(-mdot_f, 0) * phi[F_idx];
		}
		else
		{
			
			//Add contribution to owner cell
			PetscScalar val_aCC = fmax(mdot_f, 0.0);
			PetscScalar val_aCF = -fmax(-mdot_f, 0.0);

			PetscCall(MatSetValues(*A, 1, &Csol_idx, 1, &Csol_idx, &val_aCC, ADD_VALUES));
			PetscCall(MatSetValues(*A, 1, &Csol_idx, 1, &Fsol_idx, &val_aCF, ADD_VALUES));

			/* A[IDX(Csol_idx, Csol_idx, NSOLCELLS)] += fmax(mdot_f, 0.0);
			A[IDX(Csol_idx, Fsol_idx, NSOLCELLS)] += -fmax(-mdot_f, 0.0); */

			//Assign contributions to neighbor cell
			PetscScalar val_aFF = fmax(-mdot_f, 0.0); 
			PetscScalar val_aFC = -fmax(mdot_f, 0.0);

			PetscCall(MatSetValues(*A, 1, &Fsol_idx, 1, &Fsol_idx, &val_aFF, ADD_VALUES));
			PetscCall(MatSetValues(*A, 1, &Fsol_idx, 1, &Csol_idx, &val_aFC, ADD_VALUES));

			/* A[IDX(Fsol_idx, Fsol_idx, NSOLCELLS)] += fmax(-mdot_f, 0.0);
			A[IDX(Fsol_idx, Csol_idx, NSOLCELLS)] += -fmax(mdot_f, 0.0);
 */
			//Assign Source term/RHS contributions (replace with
			// function that changes for a selected scheme)
			/* b[Csol_idx] += 0.0;
			b[Fsol_idx] += 0.0; */

			// Calculate phi at the face based on higher order scheme
			double phi_f_U = phi2face(phi[C_idx], phi[F_idx], mdot_f, grad, cell_C, cell_F, f, UPWIND); // upwind scheme face value
			double phi_f_HO = phi2face(phi[C_idx], phi[F_idx], mdot_f, grad, cell_C, cell_F, f, ADVECTION_SCHEME); //high order face value

			// Add high order scheme contribution to cell RHS
			PetscScalar val_bC = -mdot_f * (phi_f_HO - phi_f_U);
			PetscScalar val_bF = mdot_f * (phi_f_HO - phi_f_U);

			PetscCall(VecSetValues(*b, 1, &Csol_idx, &val_bC, ADD_VALUES));
			PetscCall(VecSetValues(*b, 1, &Csol_idx, &val_bF, ADD_VALUES));

/* 			b[Csol_idx] += -mdot_f * (phi_f_HO - phi_f_U);
			b[Fsol_idx] += mdot_f * (phi_f_HO - phi_f_U); */
		}

	}
	return 0;
}

/* ---------------------------- grad2face routine ---------------------------- 
*  This routine computes the gradient at a face centroid (eg gradient vector) 
*  based on the gradients at the owner and neighbor cell centroids. It uses 
*  interpolation based on the normal distances to the face and a correction
*  based on forcing the vector component along the CF vector to be equal to a
*  finite difference phiF - phiC / dCF.
 -------------------------------------------------------------------------- */
int grad2face(double* grad_face, double* grad_C, double* grad_F, double* rCF, double dCF, double phi_C, double phi_F, cell* cell_C, cell* cell_F, face* f)
{
	// Average gradient weights
	//double weight_C = 1. / fmax(cell_C->volume,1.e-10); // floor to prevent divide by zero
	//double weight_F = 1. / fmax(cell_F->volume,1.e-10);
	//double denom = weight_C + weight_F;
	
	// Compute Face interplation weight
	double dCf[3] = {
		f->xc - cell_C->xc,
		f->yc - cell_C->yc,
		f->zc - cell_C->zc
	}; // position vector from cell C centroid to face f centroid

	double dfF[3] = {
		cell_F->xc - f->xc,
		cell_F->yc - f->yc,
		cell_F->zc - f->zc,
	}; // position vector from face f centroid to cell F centroid

	double S[3] = {
		f->Sx,
		f->Sy,
		f->Sz
	}; // Face area vector

	double S_mag = mag(S); //face area magnitude

	double ef[3] = {
		S[0]/S_mag,
		S[1]/S_mag,
		S[2]/S_mag
	}; //Face normal vector 

	// Face interpolation factor based on normal distances
	double weight_C = dot(dCf, ef) / (dot(dCf, ef) + dot(dfF, ef));
	double weight_F = 1 - weight_C;

	// Compute average gradient at face
	double grad_face_avg[3] = {
		(grad_C[0]*weight_C + grad_F[0]*weight_F),
		(grad_C[1]*weight_C + grad_F[1]*weight_F),
		(grad_C[2]*weight_C + grad_F[2]*weight_F)
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

// Function to interpolate a scalar known at cell centroids C and F to the intersection f'
int scal2face(double* scal_face, face* f, cell* cell_C, cell* cell_F, double scal_C, double scal_F, double* rCF, double dCF)
{
	// Compute Face interplation weight
	double dCf[3] = {
		f->xc - cell_C->xc,
		f->yc - cell_C->yc,
		f->zc - cell_C->zc
	}; // position vector from cell C centroid to face f centroid

	double dfF[3] = {
		cell_F->xc - f->xc,
		cell_F->yc - f->yc,
		cell_F->zc - f->zc,
	}; // position vector from face f centroid to cell F centroid

	double S[3] = {
		f->Sx,
		f->Sy,
		f->Sz
	}; // Face area vector

	double S_mag = mag(S); //face area magnitude

	double ef[3] = {
		S[0]/S_mag,
		S[1]/S_mag,
		S[2]/S_mag
	}; //Face normal vector 

	// Face interpolation factor based on normal distances
	double weight_C = dot(dCf, ef) / (dot(dCf, ef) + dot(dfF, ef));
	double weight_F = 1 - weight_C;

	// Compute scalar at face f' based on interpolation
	*scal_face = scal_C * weight_C + scal_F * weight_F;

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
			phi[phi_face_idx] = b->data.phi_b(b,f,0.0); // Set phi at owner cell to boundary value

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

			// Reduces to phi_b = phi_C for q_b = 0 (zero flux/outlet condition)
			double q_b = b->data.q_b(b, f, 0.0);

			if (GAMMA != 0.0)
			{
				phi[phi_face_idx] = (GAMMA * gDiff * phi[phi_owner_idx] - q_b)
					/ (GAMMA * gDiff); // Eq 8.42 from textbook
			}
			else // for advection only case
			{
				phi[phi_face_idx] = phi[phi_owner_idx];
			}
			
			break;
		}
		case Robin: {
			// For Robin, we need to calculate the equivalent boundary value based on the given phi_b, q_b, and h_infty
			// Ignore cross diffusion term for initialization? maybe. Dont have gradient at face yet

			double h_inf = b->data.robin.h_inf(b, f, 0.0);
			double phi_inf = b->data.robin.phi_inf(b, f, 0.0);;

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

// computes maximum % change for solution cells
int maxChng(double* phi, double* phi_old, int* NCELLS, int* NDEGEN_CELLS, double* epsilon)
{
	double max_change = 0.0;

	for (int i = *NDEGEN_CELLS; i < *NCELLS; i++)
	{
		double change = fabs((phi[i] - phi_old[i]) / fmax(phi[i],1e-15));
		if (change > max_change)
		{
			max_change = change;
		}
	}

	*epsilon = max_change;
	return 0;
}

double phi2face(double phi_owner, double phi_neighbor, double mdot_f,
	double* grad, cell* owner, cell* neighbor, face* f, advectionScheme scheme)
{
	double phi_f_tilde = 0.0; //normalized phi value at face
	double phi_f = 0.0; //phi value at face

	// Upwind, Downwind, and far upwind cells and phi values (far upwind is not real cell)
	double phi_C;
	double phi_D;
	double phi_U;
	cell* C;
	cell* D;

	//Determine upwind (phi_C), downwind (phi_D), and far upwind (phi_U) cell
	if (mdot_f == 0 || phi_owner == phi_neighbor)
	{
		// if there is no mass flux or if phi_C = phi_D, then 
		// there should be no change for phi at the face
		phi_f = phi_owner;
		return phi_f;
	}

	if (mdot_f > 0)
	{
		// if mass flow is positive (flow going out from owner, owner is upstream)
		phi_C = phi_owner;
		phi_D = phi_neighbor;
		C = owner;
		D = neighbor;
	}
	else
	{
		// otherwise mass flow is negative (flow going from neighbor to owner
		// neighbor is upstream)
		phi_C = phi_neighbor;
		phi_D = phi_owner;
		C = neighbor;
		D = owner;
	}

	// Calculate dCD vector
	double dCD[3] = {
		D->xc - C->xc,
		D->yc - C->yc,
		D->zc - C->zc };

	//get gradient at upwind cell
	double grad_C[3] = {
		grad[3 * C->id],
		grad[3 * C->id + 1],
		grad[3 * C->id + 2]
	};

	// Extrapolate to get far upwind value
	phi_U = phi_D - 2 * dot(grad_C, dCD);

	// Compute normalized phi_C
	double phi_C_tilde = (phi_C - phi_U) / (phi_D - phi_U);

	// Determine phi_f_tilde based on advection scheme
	switch (scheme) {
	case UPWIND: 
		phi_f = phi_C; // For upwind, we can just return phi_C without needing to compute phi_f_tilde
		return phi_f;
		//phi_f_tilde = phi_C_tilde;
		break;
	case CD:
		phi_f_tilde = 0.5 * (1 + phi_C_tilde);
		break;
	case QUICK:
		phi_f_tilde = 3. / 8. + 3. / 4. * phi_C_tilde;
		break;
	case SMART:
		if (phi_C_tilde >= 0.0 && phi_C_tilde <= 5. / 6.)
		{
			phi_f_tilde = 3. / 8. + 3. / 4. * phi_C_tilde;
		}
		else if (5. / 6. < phi_C_tilde && phi_C_tilde <= 1.)
		{
			phi_f_tilde = 1.0; 
		}
		else
		{
			phi_f_tilde = phi_C_tilde;
		}
		break;
	case BOUNDED_CD:
		if (0. <= phi_C_tilde && phi_C_tilde <= 1.0)
		{
			phi_f_tilde = 0.5 * phi_C_tilde + 0.5;
		}
		else
		{
			phi_f_tilde = phi_C_tilde;
		}
		break;
	case SOU:
		phi_f_tilde = 3. / 2. * phi_C_tilde;
		break;
	default:
		printf(stderr, "Error: Unknown advection scheme: %d", scheme);
		return 1;
	}

	// Compute phi_f based on phi_f_tilde 
	phi_f = phi_f_tilde * (phi_D - phi_U) + phi_U;
	return phi_f;
}

int calc_Residual(Mat* A, Vec* b, double* phi, cell* cells, face* faces, int* NCELLS, int* NDEGEN_CELLS, int* NFACES, double* scaled_residual)
{
	// Residual of equation Ax=b, should go to zero as solution converges. Can be used to check for convergence and also for debugging to make sure residual is decreasing after each iteration. Note that this is not the same as epsilon which is the maximum % change in phi values between iterations, but they should be correlated.
	double residual = 0.0;

	double scaling_factor = 0.0; //scale residual globally

	int NSOLCELLS = (*NCELLS) - (*NDEGEN_CELLS); // Number of solution cells (non-boundary/degenerate cells)
	
	Vec phi_vec;
	Vec residual_vec;
	PetscCall(VecCreate(PETSC_COMM_WORLD, &phi_vec));
	PetscCall(VecSetSizes(phi_vec, PETSC_DECIDE, NSOLCELLS));
	PetscCall(VecSetFromOptions(phi_vec));
	PetscCall(VecDuplicate(phi_vec, &residual_vec)); // Create vector fore residual
	
	// Set values of phi vector
	PetscScalar *phi_array; 
	PetscCall(VecGetArray(phi_vec, &phi_array)); //set phi_vec to point to phi_array

	for (int i = 0; i < NSOLCELLS; i++)
	{
		phi_array[i] = phi[i + *NDEGEN_CELLS];
	}

	PetscCall(VecRestoreArray(phi_vec, &phi_array));

	// Compute A*phi and store in residual_vec
	PetscCall(MatMult(*A, phi_vec, residual_vec));
	// Subtract b from residual_vec to get A*phi - b
	PetscCall(VecAXPY(residual_vec, -1.0, *b));
	PetscCall(VecAbs(residual_vec)); // Take absolute value of residual vector
	PetscScalar res_sum;
	PetscCall(VecSum(residual_vec, &res_sum)); // Sum of absolute values

	// Compute Scaling factor as sum of absolute values of aC*phi for all cells, where aC is the diagonal coefficient for each cell in A. This will help to prevent issues with very small or large residual values when scaling.
	// Create vector to store diagonal of A matrix
	Vec diagA;
	PetscCall(VecCreate(PETSC_COMM_WORLD, &diagA));
	PetscCall(VecSetSizes(diagA, PETSC_DECIDE, NSOLCELLS));
	PetscCall(VecSetFromOptions(diagA));

	// Get diagonal of A matrix
	PetscCall(MatGetDiagonal(*A, diagA));

	// Multiply diagonal of A with phi
	PetscCall(VecAXPY(diagA, 1.0, phi_vec));

	// Sum absolute values to get scaling factor
	PetscScalar scaling_sum;
	PetscCall(VecAbs(diagA)); // Take absolute value of diagA vector
	PetscCall(VecSum(diagA, &scaling_sum)); // Sum of absolute values

	/* 
	for (int i=0; i < NSOLCELLS; i++)
	{
		cell* C = &cells[i + (*NDEGEN_CELLS)]; //current cell, adjust index to account for degenerate cells at beginning of cells array

		double aC = A[IDX(i, i, NSOLCELLS)]; // Diagonal coefficient for cell C
		double aFphiF_sum = 0.0;

		// loop over faces to sum a_i*phi_i for the cell C
		for (int j = 0; j < C->num_faces; j++)
		{
			face* f = &faces[C->face_ids[j]];

			if (f->owner == C->id)
			{
				if (f->boundary_face)
				{
					// For boundary faces, the neighbor cell is not included in the solution and should not be included in the residual calculation since its value is determined by the boundary condition, not by solving the linear system. So we can skip this face for the residual calculation.
					continue;
				}
				double aF = A[IDX(i, f->neighbor - (*NDEGEN_CELLS), NSOLCELLS)]; // Off-diagonal coefficient for neighbor cell F
				aFphiF_sum += aF * phi[f->neighbor];
			}
			else if (f->neighbor == C->id)
			{
				double aF = A[IDX(i, f->owner - (*NDEGEN_CELLS), NSOLCELLS)]; // Off-diagonal coefficient for neighbor cell F
				aFphiF_sum += aF * phi[f->owner];
			}
		}


		double res_i = aFphiF_sum + aC * phi[C->id] - b[i]; // Residual for cell C
		residual += fabs(res_i);
		scaling_factor += fabs(aC * phi[C->id]);
	} */

	*scaled_residual = res_sum / fmax(scaling_sum, 1e-10); // Scale residual to prevent issues with very small or large values

	return 0;
}

/* -------------------------------------------------------------------------- */
/* Transient Routines */
/* -------------------------------------------------------------------------- */
int calc_time_step(cell* cells, Mat* A, int* NCELLS, int* NDEGEN_CELLS, double* t, double* min_dt)
{
	// If implicit, just return user specified timestep
	if (!EXPLICIT)
	{
		*min_dt = DT;
		if (T_FINAL - (*t) < *min_dt)
		{
			*min_dt = T_FINAL - (*t); // Adjust time step to not exceed final time
		}
		return 0;
	}

	*min_dt = INFINITY; 
	int NSOLCELLS = (*NCELLS) - (*NDEGEN_CELLS);
	
	// Create vector to store diagonal of A matrix
	Vec diagA;
	PetscCall(VecCreate(PETSC_COMM_WORLD, &diagA));
	PetscCall(VecSetSizes(diagA, PETSC_DECIDE, NSOLCELLS));
	PetscCall(VecSetFromOptions(diagA));

	PetscCall(MatGetDiagonal(*A, diagA));

	for (int i = 0; i < NSOLCELLS; i++)
	{
		cell* C = &cells[i + (*NDEGEN_CELLS)];
		PetscScalar aC;
		PetscCall(VecGetValues(diagA, 1, &i, &aC));
		
		//double aC = A[IDX(i, i, NSOLCELLS)]; // Diagonal coefficient for cell C

		
		if (aC > 0)
		{
			double dt_C = CFL* RHO * C->volume / aC; // Time step based on diffusion stability criterion for cell C
			if (dt_C < *min_dt)
			{
				(*min_dt) = dt_C;
			}
		}
	}

	if (*min_dt == INFINITY)
	{
		// If min_dt is still infinity, it means that all aC values were zero or negative
		fprintf(stderr, "Warning: All diagonal coefficients are zero or negative. Setting time step to default value of 1e-3.\n");
		return -1; 
	}

	if (T_FINAL - (*t) < *min_dt)
	{
		*min_dt = T_FINAL - (*t); // Adjust time step to not exceed final time
	}

	return 0;
}

int build_transient(Mat* A, Vec* b, double* phi, cell* cells, int* NCELLS, int* NDEGEN_CELLS, double dt)
{
	int NSOLCELLS = (*NCELLS) - (*NDEGEN_CELLS);
	for (int i = 0; i < NSOLCELLS; i++)
	{
		cell* C = &cells[i + (*NDEGEN_CELLS)];
		double fluxC = RHO * C->volume / dt;
		double fluxC_old = -RHO * C->volume / dt;

		PetscScalar val_aCC = fluxC;
		PetscCall(MatSetValues(*A, 1, &i, 1, &i, &val_aCC, ADD_VALUES));
		//A[IDX(i, i, NSOLCELLS)] += fluxC; // Add contribution to diagonal coefficient for cell C
		
		PetscScalar val_bC = -fluxC_old * phi[C->id];
		PetscCall(VecSetValues(*b, 1, &i, &val_bC, ADD_VALUES));
		//b[i] += -fluxC_old * phi[C->id]; // Add contribution to source term for cell C

	}
	return 0;
}

int explicit_update(Mat* A, Vec* b, double* phi, cell* cells, face* faces, double* dt, int* NCELLS, int* NDEGEN_CELLS)
{
	int NSOLCELLS = (*NCELLS) - (*NDEGEN_CELLS);

	double* phi_new = malloc((NEQNS * *NCELLS) * sizeof(double));
	if (!phi_new)
	{
		// Print error message to stderr stream and exit
		fprintf(stderr, "Error: Memory allocation failed for phi_old array.\n");
		return 1; // Exit with error code
	}
	phi_new = memcpy(phi_new, phi, (NEQNS * *NCELLS) * sizeof(double));

	// Create vector to store diagonal of A matrix
	Vec phi_vec;
	Vec L_vec; 
	Vec phi_new_vec;
	PetscCall(VecCreate(PETSC_COMM_WORLD, &phi_vec));
	PetscCall(VecSetSizes(phi_vec, PETSC_DECIDE, NSOLCELLS));
	PetscCall(VecSetFromOptions(phi_vec));
	PetscCall(VecDuplicate(phi_vec, &L_vec)); // Create solution vector xp with same size as bp
	PetscCall(VecDuplicate(phi_vec, &phi_new_vec)); // Create solution vector xp with same size as bp

	// Set values of phi vector
	PetscScalar *phi_array; 
	PetscCall(VecGetArray(phi_vec, &phi_array)); //set phi_vec to point to phi_array

	for (int i = 0; i < NSOLCELLS; i++)
	{
		phi_array[i] = phi[i + *NDEGEN_CELLS];
	}

	PetscCall(VecRestoreArray(phi_vec, &phi_array));

	// Calculate L, spatial operator term
	PetscCall(MatMult(*A, phi_vec, L_vec));
	PetscCall(L_vec, -1.0, *b);

	// Calculate phi new 
	PetscCall(VecWAXPY(phi_new_vec, -(*dt)/(RHO), L_vec, phi_vec)); // phi_new_vec = -dt/(RHO*C->volume)*L_vec + phi_vec
	
	for (int i = 0; i < NSOLCELLS; i++)
	{
		cell* C = &cells[i + (*NDEGEN_CELLS)]; //current cell, adjust index to account for degenerate cells at beginning of cells array

		PetscScalar phi_new_i;
		PetscScalar L_i;
		PetscCall(VecGetValues(L_vec, 1, &i, &L_i));
		PetscCall(VecGetValues(phi_new_vec, 1, &i, &phi_new_i));
		phi_new[C->id] = -(*dt) * L_i / (RHO * C->volume) + phi[C->id]; // Explicit Euler

		/* PetscScalar aC;
		PetscCall(VecGetValues(diagA, 1, &i, &aC));
		//double aC = A[IDX(i, i, NSOLCELLS)]; // Diagonal coefficient for cell C
		double aFphiF_sum = 0.0;

		// loop over faces to sum a_i*phi_i for the cell C
		// Sum over the row
		// Loop over A and calculate L = sum of all non transient terms (spatial operator)
		for (int j = 0; j < NSOLCELLS; j++)
		{
			if (j == i)
			{
				continue; // Skip diagonal term
			}
			double aF = A[IDX(i, j, NSOLCELLS)]; // Off-diagonal coefficient for neighbor cell j
			aFphiF_sum += aF * phi[j + (*NDEGEN_CELLS)]; // Add contribution from neighbor cell j to sum
		}

		double L = aC * phi[C->id] + aFphiF_sum - b[i]; //spatial operator (all discretized non transient terms) */

		//phi_new[C->id] = phi[C->id] - *dt * L / (RHO * C->volume); // Explicit Euler update for phi at cell C
	}

	// Update phi array with new values
	memcpy(phi, phi_new, (NEQNS * *NCELLS) * sizeof(double));
	free(phi_new);
	return 0;
}

/* -------------------------------------------------------------------------- */
/*                       Known Scalar Operator Routines                       */
/* -------------------------------------------------------------------------- */
/* --------------------------- Numerical Gradient --------------------------- */
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
		double dxk = F->xc - C->xc;
		double dyk = F->yc - C->yc;
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

/* -------------------------- Numerical Divergence -------------------------- 
* This routine calculates the numerical divergence of a given scalar array phi
* and adds it to the RHS b vector for that equation. The divergence array is not stored.
* Scalar components of vectors should be contiguous in memory, so double* phi should be 
* ordered as [phi_x for all cells, phi_y for all cells, phi_z for all cells] and indexed
* using IDX macro. 
* This will be called for the continuity equation's RHS b vector, but it will use the
* guessed momentum eqation's phi values and gradients to calcate the divergence of velocity
 -------------------------------------------------------------------------- */
int build_div(Vec* b, double** phi, double** grad, int eq_ids[static ND], node* nodes, cell* cells, face* faces, int* NCELLS, int* NDEGEN_CELLS, int* NFACES)
{
	int NSOLCELLS = (*NCELLS) - (*NDEGEN_CELLS);
	/* Vec div_phi; 

	// Create vector to store divergence of phi vector
	PetscCall(VecCreate(PETSC_COMM_WORLD, &div_phi));
	PetscCall(VecSetSizes(div_phi, PETSC_DECIDE, NSOLCELLS));
	PetscCall(VecSetFromOptions(div_phi));

	// Set all entries to zero 
	PetscCall(VecZeroEntries(div_phi));	 */

	// Loop over all faces and add contributions to the b vector for owner cell C and neighbor cell F
	for (int i = 0; i < *NFACES; i++)
	{
		
		face* f = &faces[i]; //current face 

		cell* C = &cells[f->owner]; //owner cell
		cell* F = &cells[f->neighbor]; //neighbor cell

		// Calculate position vectors/ Geometry stuff
		double rCF[3] = {
			F->xc - C->xc,
			F->yc - C->yc,
			F->zc - C->zc
		}; // position vector from cell C centroid to cell F centroid

		double dCF = mag(rCF); // distance between cell centroids

		// Calculate rf' position vector
		double Sf[3] = {f->Sx, f->Sy, f->Sz}; //face area vector
		double Sf_mag = mag(Sf); // magnitude of face area vector
		double ef[3] = {Sf[0]/Sf_mag, Sf[1]/Sf_mag, Sf[2]/Sf_mag}; // Normal vector of face

		double rCf[3] = {
			f->xc - C->xc,
			f->yc - C->yc,
			f->zc - C->zc
		}; // position vector from cell C centroid to face f centroid

		double g_f = dot(rCf, ef)/dot(rCF, ef); // interpolation factor for gradient at face based on normal distances

		// vector to the normal intersection of the face from cell C centroid
		double rf_prime[3] = {
			(1 - g_f) * C->xc + g_f * F->xc,
			(1 - g_f) * C->yc + g_f * F->yc,
			(1 - g_f) * C->zc + g_f * F->zc
		};

		double r_f_prime_f[3] = {
			f->xc - rf_prime[0],
			f->yc - rf_prime[1],
			f->zc - rf_prime[2]
		}; // vector from face intersection to face centroid

		// Do not need a boundary face if statement since the values of the scalar at the boundary should be stored in phi
		double vec_C[3] = {0.0, 0.0, 0.0};
		double vec_F[3] = {0.0, 0.0, 0.0};
		double vec_f_prime[3] = {0.0, 0.0, 0.0}; 
		double vec_f[3] = {0.0, 0.0, 0.0};

		// velocity vector gradient
		double jacobian_C[ND][3]; // jacobian for cell c, gradient of velocity vector
		double jacobian_F[ND][3]; // jacobian for cell F, gradient of velocity vector
		double jacobian_f_prime[ND][3]; // jacobian for face f, gradient of velocity vector
		
		// Assemble vectors. loop over components
		for (int dims = 0; dims < ND; dims++)
		{
			vec_C[dims] = phi[eq_ids[dims]][C->id];
			vec_F[dims] = phi[eq_ids[dims]][F->id];

			// Interpolate the vector to face intersection
			scal2face(&vec_f_prime[dims], f, C, F, vec_C[dims], vec_F[dims], rCF, dCF);

			// interpolate gradient of each vector component
			// get the gradient component (all ways 3) at the centroid
			jacobian_C[dims][0] = grad[eq_ids[dims]][IDX(0, C->id, 3)];
			jacobian_C[dims][1] = grad[eq_ids[dims]][IDX(1, C->id, 3)];
			jacobian_C[dims][2] = grad[eq_ids[dims]][IDX(2, C->id, 3)];

			jacobian_F[dims][0] = grad[eq_ids[dims]][IDX(0, F->id, 3)];
			jacobian_F[dims][1] = grad[eq_ids[dims]][IDX(1, F->id, 3)];
			jacobian_F[dims][2] = grad[eq_ids[dims]][IDX(2, F->id, 3)];

			// interpolate gradients to face
			grad2face(&jacobian_f_prime[dims], jacobian_C[dims], jacobian_F[dims], rCF, dCF, vec_C[dims], vec_F[dims], C, F, f);

			// Correct face values for skewness
			vec_f[dims] = vec_f_prime[dims] + dot(jacobian_f_prime[dims], r_f_prime_f); 
		}

		// Calculate and add RHS vector contributions
		PetscScalar val_bC = dot(vec_f, Sf); // contribution to owner cell C
		 PetscCall(VecSetValues(*b, 1, &C->id, &val_bC, ADD_VALUES));
		//b[C->id] += dot(vec_f, Sf); // contribution to owner

		PetscScalar val_bF = -dot(vec_f, Sf); // contribution to neighbor cell F (negative since normal vector points outward from owner cell)
		 PetscCall(VecSetValues(*b, 1, &F->id, &val_bF, ADD_VALUES));
		//b[F->id] += -dot(vec_f, Sf); // contribution to neighbor cell F
	}	

	return 0;
}

/* --------------------------- Numerical Laplacian -------------------------- 
* Numerical laplacian operator. Operates on a known scalar field phi
* and adds contributions to the RHS b vector for that equation. The laplacian
* array is not stored. Not technically laplacian becasue it also includes a 
* diffusion coefficient for adding to the RHS. For a pure laplacian operator,
* the diffusion coefficient should be set to 1.0.
 -------------------------------------------------------------------------- */
int build_lap(Vec* b, double* phi, node* nodes, cell* cells, face* faces, boundary* boundaries, int* NCELLS, int* NDEGEN_CELLS, int* NFACES)
{

}
