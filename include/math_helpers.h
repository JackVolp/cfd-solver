#ifndef MATH_HELPERS_H
#define MATH_HELPERS_H


#define dot(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
//magnitude macro, a lot easier to use than the function 
#define mag(a) (sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])) 

void magnitude(double* v, double* result);

void cross_prod(double* a, double* b, double* result);

void solve_2x2_system(double A11,
	double A12,
	double A21,
	double A22,
	double b1,
	double b2,
	double* x1,
	double* x2);


#endif // !MATH_HELPERS_H
