
#ifndef SETUP_H
#define SETUP_H

#include <stdbool.h>
#include "grid.h"
//#include "solver.h"


// Solver parameters
#define STOP_COND 1e-15
#define RPRT_INTERVAL 1
#define MAX_ITER 10

// Possible advection Schemes user can choose 
typedef enum advectionScheme {
	UPWIND = 0,
	CD = 1,
	QUICK = 2,
	SMART = 3,
	BOUNDED_CD = 4,
	SOU = 5
} advectionScheme;

extern advectionScheme ADVECTION_SCHEME;

// Problem Setup
#define transient false
#define CFL 0.5 // CFL number, only used when transient
#define T_FINAL 1.0 // Final time for transient simulation
#define NEQNS 1 // Number of transport equations solved


// Physical constants
//Thermal Diffusivity of aluminum at room temperature (m^2/s)
//#define GAMMA 1
#define GAMMA 0 //set diffusion to zero for pure advection
#define RHO 1 //Density

// Source Term
//#define Q_C(x,y,z) (10.*x + 5.)
#define Q_C(x,y,z) (0)

// Velocity Field
#define XVEL 1.0 //x-velocity
#define YVEL 1.0 //y-velocity

// Boundary profiles
double inlet_profile(const boundary* b, const face* f, double t);
double phi0_boundary(const boundary* b, const face* f, double t);
double zero_flux(const boundary* b, const face* f, double t);

#endif // !SETUP_H