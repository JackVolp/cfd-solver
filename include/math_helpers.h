#ifndef MATH_HELPERS_H
#define MATH_HELPERS_H

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
