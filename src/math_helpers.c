#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>


void magnitude(double* v, double* result)
{
	*result = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void cross_prod(double* a, double* b, double* result)
{
	// indexing is equal to a[i] == *(a + i)
	result[0] = a[1] * b[2] - a[2] * b[1];
	result[1] = a[2] * b[0] - a[0] * b[2];
	result[2] = a[0] * b[1] - a[1] * b[0];
}

void solve_2x2_system(double A11, double A12, double A21, double A22,
	double b1, double b2,
	double* x1, double* x2)
{
	double det = A11 * A22 - A12 * A21;
	if (det == 0)
	{
		fprintf(stderr, "Error: Coefficient matrix is singular, cannot solve system.\n");
		exit(1);
	}

	*x1 = (A22 * b1 - A12 * b2) / det;
	*x2 = (A11 * b2 - A21 * b1) / det;
}