#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cfd.h"

/* Indexing macros (row major)*/
//#define IDX(i,j,nx) ((j)*(nx) + (i))
#define IDX(i,j,eq,nx,ny) ((j)*(nx) + (i) + (eq)*(nx*ny))

int main(void)
{

	// Create Read grid function to read the points and then create arrays 
	// Grid Setup
	const int NX = 100; //Number Divisions x-direction
	const int NY = 100; //Number divisions y-direction
	const int NZ = 1;   //Number divisions z-direction 

	//Domain setup
	const double x0 = 0.0, y0 = 0.0, z0 = 0.0; // Origin Location
	const double x1 = 1.0, y1 = 1.0, z1 = 1.0; // Domain side lengths

	//Grid Spacing (only 2D array for now)
	const double dx = (x1 - x0) / NX; // Grid spacing x-direciton
	const double dy = (y1 - y0) / NY; // Grid spacing y-direction
	const double dz = (z1 - z0) / NZ; // Grid spacing z-direction

	const int NCELLS = NX * NY * NZ;
	const int NPOINTS = (NX + 1) * (NY + 1) * (NZ + 1);

	// Allocate Arrays
	const int NEQNS = 1; // Number of equations solved

	// Conservative variable fields [rho, rho*u, rho*v]
	double* phi = malloc((NCELLS * NEQNS)* sizeof(double)); // Solution variable

	
	// Check for correct memory allocation
	if (phi == NULL)
	{
		// Print error message to stderr stream and exit
		fprintf(stderr, "Error: Memory allocation failed for phi array.\n");
		return 1; // Exit with error code
	}

	// Nodes array
	node* nodes = malloc(NPOINTS * sizeof(node)); // Node coordinates array

	if (nodes == NULL)
	{
		fprintf(stderr, "Error: Memory allocation failed for nodes array.\n");
		return 1; // Exit with error code
	}

	// Calculate Node Coordinates
	for (int i = 0; i < NPOINTS; i++)
	{

	}
	// Loop over grid and initialize conservative ariables
	for (int eq = 0; eq < NEQNS; eq++) //loop over conservative fields
	{
		for (int j = 0; j < NY; j++) // loop over y
		{
			for (int i = 0; i < NX; i++) //loop over x
			{
				//Global index
				int idx = IDX(i, j, eq, NX, NY);

				// Cell Centroid coordinates
				double x = x0 + dx / 2 + i * dx; //x cell centroid points
				double y = y0 + dy / 2 + j * dy; //y cell centroid points 
				
				// Allocate Fields
				if (eq == 0) { phi[idx] = 1; } else //Density
				if (eq == 1) { phi[idx] = 0; } else //x-momentum
				if (eq == 2) { phi[idx] = 0; } else //y-momentum
				{ phi[idx] = 0; }                   //other unspecified 


			}
		}
	}
	

	/*----- Write legacy .vtk file-------*/
	// Create file 
	FILE* fp = fopen("uniform_2d.vtk", "w"); //open file in write mode
	if (!fp)
	{
		perror("Error writing data file\n");
		return 1;
	}

	// Header section
	fprintf(fp, "# vtk DataFile Version 3.0\n");
	fprintf(fp, "Jack's first vtk datafile\n");
	fprintf(fp, "ASCII\n");

	// Data definition
	fprintf(fp, "DATASET STRUCTURED_POINTS\n");

	// Geometry
	fprintf(fp,"DIMENSIONS %d %d %d\n", NX+1, NY+1, NZ+1); //Number of grid points (not cells)
	fprintf(fp,"ORIGIN %f %f %f\n", x0, y0, z0); 
	fprintf(fp,"SPACING %f %f %f\n", dx, dy, dz);

	// Data
	fprintf(fp, "CELL_DATA %d\n", NCELLS);
	fprintf(fp, "SCALARS phi(0) float 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");

	// Save Density Scalar Data
	for (int j = 0; j < NY; j++)
	{
		for (int i = 0; i < NX; i++)
		{
			int idx = IDX(i, j, 0, NX, NY);
			fprintf(fp, "%f\n", (float)phi[idx]);
		}
	}
	
	/*
	// Save Momentum Vector Data
	fprintf(fp, "VECTORS momentum float\n");
	for (int j = 0; j < NY; j++)
	{
		for (int i = 0; i < NX; i++)
		{
			int idx_mx = IDX(i, j, NX, NY, 1);
			int idx_my = IDX(i, j, NX, NY, 2);

			fprintf(fp, "%f %f %f\n", (float)phi[idx_mx], (float)phi[idx_my], 0.0f);
		}
	}
	*/

	fclose(fp);

	// Free Memory
	free(phi);
	free(nodes);

	printf("To C or not to C: that is the question. \n");
	return 0;
}