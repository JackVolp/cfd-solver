
#ifndef SETUP_H
#define SETUP_H

#include <stdbool.h>
#include "config.h"
#include "grid.h"
//#include "solver.h"

/* -------------------------------------------------------------------------- */
/*                                     I/O                                    */
/* -------------------------------------------------------------------------- */
/* -------------------------- Grid input file name -------------------------- */
extern const char *filename;

/* ---------------------------- Output file name ---------------------------- */
extern const char *out_fname;
/* -------------------------------------------------------------------------- */

// Solver parameters
#define STOP_COND 1e-15
#define RPRT_INTERVAL 10
#define MAX_ITER 100000
#define C_EPSILON 0.5 // stabilization constant for Brezzi Pitkaranta stabilization O(1)
// Problem Setup
// For a transient problem, I need:
	// T_FINAL 
	// if explicit:
	// CFL condition to calculate time step (13.27)
	// if implicit:
	// time step size 

#define TRANSIENT false
#define EXPLICIT false
#define CFL 0.2 // CFL number, only used when transient and explicit are both true
#define DT 0.5 // time step size, only used when transient is true but explicit is false (inplicit stepping tstep size)
#define T_FINAL 7.5 // Final time for transient simulation only used when transient
#define SAVE_INTERVAL 0.5 // Time interval for saving output files, only used when transient

#define SIMPLE true //Use simple for p-v coupling. If yes, adittional pressure array is allocated since continuity equation is solved with a correction

// Physical constants
//Thermal Diffusivity of aluminum at room temperature (m^2/s)
//#define GAMMA 1
//#define GAMMA 0 //set diffusion to zero for pure advection
#define RHO 1 //Density

//Velocity field
#define XVEL 1.0
#define YVEL 1.0

/* ---------------------------- Equation indicies --------------------------- */
#define XMOM 0
#define YMOM 1
#define PCORR 2

// Possible advection Schemes user can choose 
typedef enum advectionScheme {
	NONE = -1,
	UPWIND = 0,
	CD = 1,
	QUICK = 2,
	SMART = 3,
	BOUNDED_CD = 4,
	SOU = 5
} advectionScheme;

extern advectionScheme ADVECTION_SCHEME;

typedef double diffusionCoeff;
extern diffusionCoeff GAMMA[NEQNS];

// Boundary profiles
double no_slip_wall(const boundary* b, const face* f, double t);
/* double inlet_profile(const boundary* b, const face* f, double t);
double phi0_boundary(const boundary* b, const face* f, double t); */
double zero_flux(const boundary* b, const face* f, double t);

// Boundary Conditions
#define NBOUNDARIES 4 // Number of unique boundary conditions for the problem
extern boundaryType problem_boundary_types[NBOUNDARIES][NEQNS];
extern boundaryData problem_boundary_data[NBOUNDARIES][NEQNS];

/* ------------------------- Source Term Definition ------------------------- */
// Source term function type definition
typedef double (*SOURCE_TERM_FCN)(const cell* C, const double t); 
extern SOURCE_TERM_FCN SOURCE_TERMS[NEQNS]; // Array of source term functions for each equation
double zero_source(const cell* C, const double t);
double fx(const cell* C, const double t);
double fy(const cell* C, const double t);
double g(const cell* C, const double t);

#endif // !SETUP_H