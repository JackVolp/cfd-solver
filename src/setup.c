#include "setup.h"
#include <math.h>

/* -------------------------------------------------------------------------- */
/*                                     I/O                                    */
/* -------------------------------------------------------------------------- */
/* -------------------------- Grid input file name -------------------------- */
const char *filename = "input/hw2_64x64.vtk";

/* ---------------------------- Output file name ---------------------------- */
const char *out_fname = "output/hw2_64x64_implicit_BOUNDED_CD_out.vtk";	
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

/* -------------------------------------------------------------------------- */
/* Advection Scheme Selection */
/* -------------------------------------------------------------------------- */
//advectionScheme ADVECTION_SCHEME = UPWIND;
advectionScheme ADVECTION_SCHEME = BOUNDED_CD;