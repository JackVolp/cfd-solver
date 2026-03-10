#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "cfd.h"

/* Indexing macros (row major)*/
//#define IDX(i,j,nx) ((j)*(nx) + (i))
#define IDX(i,j,eq,nx,ny) ((j)*(nx) + (i) + (eq)*(nx*ny))

int main(void)
{		
	node* nodes;
	cell* cells; 

	int NPOINTS = 0, NCELLS = 0, CELL_LIST_SIZE = 0;

	int err = read_grid("C:\\Users\\jvolponi0552\\Documents\\GitHub\\cfd-solver\\gmsh_grid.vtk", &nodes, &cells, &NPOINTS, &NCELLS, &CELL_LIST_SIZE);

	free(nodes);
	free(cells);

	

	/*
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
	
	*/


	/*----- Write legacy .vtk file-------*/

	/*
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
	

	fclose(fp);

	// Free Memory
	free(phi);
	free(nodes);

	*/

	printf("To C or not to C: that is the question. \n");
	return 0;
	
}

int read_grid(const char* filename, node** nodes_out, cell** cells_out, int* NPOINTS, int* NCELLS, int* CELL_LIST_SIZE)
{
	// pass pointer as double pointer to allow modification of caller's pointer to nodes and cells arrays

	// Open file for reading
	FILE* fp = fopen(filename, "r");
	if (!fp)
	{
		perror("Error opening grid file\n");
		return 1; // Exit with error code
	}

	// Read number of nodes and cells
	char line[512]; //Store line from file

	bool points_section = false;
	bool cells_section = false;
	bool cell_types_section = false;

	int pidx = 0; // Node index
	int cidx = 0; // Cell index
	int ctidx = 0; // Cell type index

	// Create Null Pointers for nodes and cells, will allocate memory after reading number of points and cells
	node* nodes = NULL;
	cell* cells = NULL;

	while (fgets(line, sizeof(line), fp))
	{
		if (line[0] == '\n' || line[0] == '\r') // Skip empty lines and comments
			continue;

		// Check for section headers
		if (sscanf(line, "POINTS %d", NPOINTS) == 1)
		{
			points_section = true;
			cells_section = false;
			cell_types_section = false;

			// Allocate Nodes memory after reading number of points
			nodes = malloc((size_t)(*NPOINTS) * sizeof(node));
			if (nodes == NULL)
			{
				fprintf(stderr, "Error: Memory allocation failed for nodes array.\n");
				fclose(fp);
				return 2; // Exit with error code
			}

			pidx = 0; // Reset node index for reading node data
			continue;
		}

		if (sscanf(line, "CELLS %d %d", NCELLS, CELL_LIST_SIZE) == 2)
		{
			cells_section = true;
			points_section = false;
			cell_types_section = false;

			//Allocate Cells memory after reading number of cells
			cells = malloc((size_t)(*NCELLS) * sizeof(cell));
			if (cells == NULL)
			{
				fprintf(stderr, "Error: Memory allocation failed for cells array.\n");
				fclose(fp);
				return 2; // Exit with error code
			}

			cidx = 0; // Reset cell index for reading cell data
			continue;
		}

		if (sscanf(line, "CELL_TYPES %d", NCELLS) == 1)
		{
			cell_types_section = true;
			points_section = false;
			cells_section = false;

			ctidx = 0; // Reset cell type index for reading cell type data
			continue;
		}


		// Read node data, cell data, or cell type data based on the current section
		if (points_section)
		{
			double x, y, z;
			if (sscanf(line, "%lf %lf %lf", &x, &y, &z) == 3)
			{
				nodes[pidx].x = x;
				nodes[pidx].y = y;
				nodes[pidx].z = z;

				nodes[pidx].id = pidx;

				pidx++;
				continue;
			}
			else
			{
				fprintf(stderr, "Error reading node data\n");
				return 1;
			}
		}
		else if (cells_section)
		{
			int num_nodes;
			int node_ids[8]; // Assuming max 8 nodes per cell
			for (int i = 0; i < 8; i++)
			{
				cells[cidx].node_ids[i] = -1; // Initialize node IDs to -1
			}

			if (sscanf(line, "%d %d %d %d %d %d %d %d %d", &num_nodes, &node_ids[0], &node_ids[1], &node_ids[2], &node_ids[3], &node_ids[4], &node_ids[5], &node_ids[6], &node_ids[7]) >= 2)
			{

				// Store node IDs for the cell
				cells[cidx].num_nodes = num_nodes; // Store number of nodes in the cell
				for (int i = 0; i < num_nodes; i++)
				{
					cells[cidx].node_ids[i] = node_ids[i];

				}
				// Store cell ID
				cells[cidx].id = cidx;

				cidx++;
			}
			else
			{
				fprintf(stderr, "Error reading cell data\n");
				return 1;
			}
		}
		else if (cell_types_section)
		{
			int cell_type;
			if (sscanf(line, "%d", &cell_type) == 1)
			{
				cells[ctidx].type = cell_type;
			}
			else
			{
				fprintf(stderr, "Error reading cell type data\n");
				return 1;
			}

			ctidx++;
		}
	}

	// close file
	fclose(fp);

	//Release memory on error 
	if (nodes && pidx != *NPOINTS) { free(nodes); free(cells); return 11; }
	if (cells && cidx != *NCELLS) { free(nodes); free(cells); return 12; }

	// Set output pointers
	*nodes_out = nodes;
	*cells_out = cells;

	return 0; // Success

}