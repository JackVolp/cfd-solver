#include "setup.h"
#include <math.h>

/* -------------------------------------------------------------------------- */
/*                                     I/O                                    */
/* -------------------------------------------------------------------------- */
/* -------------------------- Grid input file name -------------------------- */
const char *filename = "input/hw3_128x128.vtk";
//const char *filename = "input/full_unstruct.vtk";

/* ---------------------------- Output file name ---------------------------- */
const char *out_fname = "output/hw3_128x128_out.vtk";	
//const char *out_fname = "output/full_unstruct_implicit_SMART_out.vtk";

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* Functions related to the problem setup. BC_PROFILES etc */
/* -------------------------------------------------------------------------- */
// Use for velocity boundary conditions
double no_slip_wall(const boundary* b, const face* f, double t)
{
	(void)t;
	(void)b;
	(void)f;

	return 0.0;
}

// Use for pressure
double zero_flux(const boundary* b, const face* f, double t)
{
	(void)t;
	(void)b;
	(void)f;

	return 0.0;
}

/* ---------------------------- Define Boundaries --------------------------- */
// Initialize boundaries (change to allocate for more complex gemoetry)
// probably move this to setup somehow also
// Define boundaries in the order they were defined in the .vtk file. Make sure interior
// has the id=9

boundaryType problem_boundary_types[NBOUNDARIES][NEQNS] = { 
	{Dirichlet, Dirichlet, Neumann}, //List of boundary 0 types for all equations 
	{Dirichlet, Dirichlet, Neumann}, //List of boundary 1 types for all equations 
	{Dirichlet, Dirichlet, Neumann}, //List of boundary 2 types for all equations 
	{Dirichlet, Dirichlet, Neumann}  //List of boundary 3 types for all equations
}; // problem_boundary_types[0] = &problem_boundary_types


/* boundaryType problem_boundary_types[NBOUNDARIES][NEQNS] = { 
	{Dirichlet}, //List of boundary 0 types for all equations 
	{Neumann},	 //List of boundary 1 types for all equations 
	{Dirichlet}	 //List of boundary 2 types for all equations 
}; // problem_boundary_types[0] = &problem_boundary_types */

boundaryData problem_boundary_data[NBOUNDARIES][NEQNS] = {
    {
        {.phi_b = no_slip_wall},  // boundary 0, equation 0
        {.phi_b = no_slip_wall},  // boundary 0, equation 1
        {.q_b   = zero_flux}      // boundary 0, equation 2
    },
    {
        {.phi_b = no_slip_wall},  // boundary 1, equation 0
        {.phi_b = no_slip_wall},  // boundary 1, equation 1
        {.q_b   = zero_flux}      // boundary 1, equation 2
    },
    {
        {.phi_b = no_slip_wall},  // boundary 2, equation 0
        {.phi_b = no_slip_wall},  // boundary 2, equation 1
        {.q_b   = zero_flux}      // boundary 2, equation 2
    },
    {
        {.phi_b = no_slip_wall},  // boundary 3, equation 0
        {.phi_b = no_slip_wall},  // boundary 3, equation 1
        {.q_b   = zero_flux}      // boundary 3, equation 2
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
// pressure correction eqn (eqn 2) gets place holder zero since its coefficient is grid dependent
diffusionCoeff GAMMA[NEQNS] = {1.0, 1.0, 0.0}; 
/* --------------------- Source Terms for each equation --------------------- */
double zero_source(const cell* C, const double t)
{
	(void)C;
	(void)t;

	return 0;
}

//x momentum eqn source term
double fx(const cell* C, const double t)
{
	(void)t;

	double pi = acos(-1.0); //pi
	double x = C->xc;
	double y = C->yc;

	double term1 = 2*pi*sin(2*pi*x);
	double term2 = GAMMA[0]*(2*pi*pi*pi*cos(2*pi*x)*sin(2*pi*y) 
		- 4*pi*pi*pi*sin(pi*x)*sin(pi*x)*sin(2*pi*y));
	
	return term1 + term2;
}

//y momentum eqn source term
double fy(const cell* C, const double t)
{
	(void)t;

	double pi = acos(-1.0); //pi
	double x = C->xc;
	double y = C->yc;

	double term1 = 2*pi*sin(2*pi*y);
	double term2 = GAMMA[1]*(4.*pi*pi*pi*sin(2.*pi*x)*sin(pi*y)*sin(pi*y) 
		- 2.*pi*pi*pi*sin(2.*pi*x)*cos(2.*pi*y));
	
	return term1 + term2;
}

// continuity eqn source term. Since alpha = epsilon, epsilon = (alpha + epsilon) / 2
double EPSILON_G = 0.0;
double g(const cell* C, const double t)
{
	(void)t;

	double pi = acos(-1.0); //pi
	double x = C->xc;
	double y = C->yc;

	//double epsilon = GAMMA[PCORR]/2.0; // remove alpha contribution from diffusion coefficient
	
	return -4.0 * EPSILON_G * pi * pi * (cos(2. * pi * x) + cos(2. * pi * y));
}

//SOURCE_TERM_FCN SOURCE_TERMS[NEQNS] = {zero_source}; // Array of source term functions for each equation, initialized to zero source for all equations
SOURCE_TERM_FCN SOURCE_TERMS[NEQNS] = {fx, fy, g}; // Array of source term functions for each equation, initialized to specific source terms for each equation