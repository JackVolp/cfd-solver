
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <stdbool.h>

// Solver parameters
#define STOP_COND 1e-15
#define RPRT_INTERVAL 1
#define MAX_ITER 100

// Problem Setup
#define transient false
#define CFL 0.5 // CFL number, only used when transient
// Physical constants
//Thermal Diffusivity of aluminum at room temperature (m^2/s)
#define GAMMA 1

// Source Term
#define Q_C(x,y,z) (10.*x + 5.)
//#define Q_C(x,y,z) (0)

#endif // !CONSTANTS_H