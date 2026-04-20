
#ifndef SETUP_H
#define SETUP_H

#include <stdbool.h>
#include <grid.h>

// Solver parameters
#define STOP_COND 1e-15
#define RPRT_INTERVAL 1
#define MAX_ITER 0

// Problem Setup
#define transient false
#define CFL 0.5 // CFL number, only used when transient
#define NEQNS 1 // Number of transport equations solved

// Physical constants
//Thermal Diffusivity of aluminum at room temperature (m^2/s)
#define GAMMA 1

// Source Term
#define Q_C(x,y,z) (10.*x + 5.)
//#define Q_C(x,y,z) (0)

// Boundary profiles
double inlet_profile(const boundary* b, const face* f, double t);
double phi0_boundary(const boundary* b, const face* f, double t);
double zero_flux(const boundary* b, const face* f, double t);

#endif // !SETUP_H