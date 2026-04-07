
#ifndef CONSTANTS_H
#define CONSTANTS_H

// Solver parameters
<<<<<<< HEAD
#define MAX_ITER 10000
#define STOP_COND 1e-15
#define RPRT_INTERVAL 1
=======
#define MAX_ITER 100
>>>>>>> 2047fb40847928e79c1cbbed2c1d77c94f81c8b0

// Physical constants
//Thermal Diffusivity of aluminum at room temperature (m^2/s)
#define GAMMA 1

// Source Term
#define Q_C(x,y,z) (10.*x + 5.)
//#define Q_C(x,y,z) (0)

#endif // !CONSTANTS_H