#include "setup.h"
#include <math.h>

/* -------------------------------------------------------------------------- */
/*                                     I/O                                    */
/* -------------------------------------------------------------------------- */
/* -------------------------- Grid input file name -------------------------- */
const char *filename = "input/hw3_64x64.vtk";
//const char *filename = "input/full_unstruct.vtk";

/* ---------------------------- Output file name ---------------------------- */
const char *out_fname = "output/hw3_64x64_out.vtk";	
//const char *out_fname = "output/full_unstruct_implicit_SMART_out.vtk";

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* Functions related to the problem setup. BC_PROFILES etc */
/* -------------------------------------------------------------------------- */
double south_profile(const boundary* b, const face* f, double t)
{
	(void)t;
	(void)b;

	return sin(f->xc);
}

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
// Define boundaries in the order they were defined in the .vtk file. Make sure interior
// has the id=9
boundaryType problem_boundary_types[NBOUNDARIES][NEQNS] = { 
	{Dirichlet}, //List of boundary 0 types for all equations 
	{Neumann},	 //List of boundary 1 types for all equations 
	{Dirichlet}	 //List of boundary 2 types for all equations 
}; // problem_boundary_types[0] = &problem_boundary_types

boundaryData problem_boundary_data[NBOUNDARIES][NEQNS] = {
	{
		{.phi_b = (*phi0_boundary)} // List of boundary 0 data for all equations
	},
	{
		{.q_b = (*zero_flux)} // List of boundary 1 data for all equations
	},
	{
		{.phi_b = (*inlet_profile)} //List of boundary 2 data for all equations
	}
};

/* boundaryData problem_boundary_data[NBOUNDARIES][NEQNS] = { 
	{ //Boundary data for equation 0
	{.phi_b = (*phi0_boundary)},
	{.q_b = (*zero_flux)},
	{.phi_b = (*inlet_profile)}
}//,
	//{ //Boundaries for equation 1 (if there were a second equation, e.g. for velocity in a coupled flow problem)
	//{.phi_b = (*phi0_boundary)},
	//{.q_b = (*zero_flux)},
	//{.phi_b = (*inlet_profile)}
	//}
}; */

/* boundaryData problem_boundary_data[NBOUNDARIES][NEQNS] {
	{}
}] */
/* -------------------------------------------------------------------------- */
/* Advection Scheme Selection */
/* -------------------------------------------------------------------------- */
//advectionScheme ADVECTION_SCHEME = UPWIND;
advectionScheme ADVECTION_SCHEME = SMART;

/***************** Diffusion Coefficients for each equation *****************/
diffusionCoeff GAMMA[NEQNS] = {0.0};

/* --------------------- Source Terms for each equation --------------------- */
double zero_source(const cell* C, const double t)
{
	(void)C;
	(void)t;

	return 0;
}
SOURCE_TERM_FCN SOURCE_TERMS[NEQNS] = {zero_source}; // Array of source term functions for each equation, initialized to zero source for all equations