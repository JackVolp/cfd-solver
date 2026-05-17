#include "setup.h"
#include <math.h>

/* -------------------------------------------------------------------------- */
/*                                     I/O                                    */
/* -------------------------------------------------------------------------- */
/* -------------------------- Grid input file name -------------------------- */
//const char *filename = "input/hw2_64x64.vtk";
const char *filename = "input/full_unstruct.vtk";

/* ---------------------------- Output file name ---------------------------- */
//const char *out_fname = "output/hw2_64x64_implicit_BOUNDED_CD_out.vtk";	
const char *out_fname = "output/full_unstruct_implicit_SMART_out.vtk";

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* Functions related to the problem setup. BC_PROFILES etc */
/* -------------------------------------------------------------------------- */
double inlet_profile(const boundary* b, const face* f, double t)
{
	(void)t;
	(void)b;

	return sin(f->yc);
}

double phi0_boundary(const boundary* b, const face* f, double t)
{
	(void)t;
	(void)b;
	(void)f;

	return 0.0;
}
double zero_flux(const boundary* b, const face* f, double t)
{
	(void)t;
	(void)b;
	(void)f;

	return 0;
}

/* ---------------------------- Define Boundaries --------------------------- */
// Initialize boundaries (change to allocate for more complex gemoetry)
// probably move this to setup somehow also
boundaryType problem_boundary_types[NBOUNDARIES] = {Dirichlet, Neumann, Dirichlet};

// boundaryType p1_boundaries[4] = { Neumann, Robin, Neumann, Dirichlet };
/*boundaryData p1_boundary_data[4] = {
	{.q_b = 0.0},
	{.robin = {.h_inf = 100.0, .phi_inf = 25.}},
	{.q_b = 0.0},
	{.phi_b = 100.}
};*/

boundaryData problem_boundary_data[NBOUNDARIES] = {
	{.phi_b = (*phi0_boundary)},
	{.q_b = (*zero_flux)},
	{.phi_b = (*inlet_profile)}};

/* -------------------------------------------------------------------------- */
/* Advection Scheme Selection */
/* -------------------------------------------------------------------------- */
//advectionScheme ADVECTION_SCHEME = UPWIND;
advectionScheme ADVECTION_SCHEME = SMART;

/***************** Diffusion Coefficients for each equation *****************/
diffusionCoeff GAMMA[NEQNS] = {0.0};