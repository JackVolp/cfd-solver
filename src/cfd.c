#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "cfd.h"

/* Indexing macros (row major)*/
//#define IDX(i,j,nx) ((j)*(nx) + (i))
//#define IDX(i,j,eq,nx,ny) ((j)*(nx) + (i) + (eq)*(nx*ny))

#define IDX(i,eq,NCELLS) ( (i) + (eq*NCELLS) )

int main(void)
{		
	// Load grid
	node* nodes;
	cell* cells; 

	int NPOINTS = 0, NCELLS = 0, CELL_LIST_SIZE = 0;

	//Home PC
	int err = read_grid("C:\\Users\\jtvol\\Documents\\ME696\\Convection-Diffusion\\out\\build\\x64-Debug\\gmsh_grid.vtk", &nodes, &cells, &NPOINTS, &NCELLS, &CELL_LIST_SIZE);
	
	//Lab PC
	//int err = read_grid("C:\\Users\\jvolponi0552\\Documents\\GitHub\\cfd-solver\\gmsh_grid.vtk", &nodes, &cells, &NPOINTS, &NCELLS, &CELL_LIST_SIZE);

	int NEQNS = 3; // Number of transport equations solved

	// Allocate conservative scalars
	double* phi = malloc((NEQNS * NCELLS) * sizeof(double));
	if (phi == NULL)
	{
		// Print error message to stderr stream and exit
		fprintf(stderr, "Error: Memory allocation failed for phi array.\n");
		return 1; // Exit with error code
	}

	// initialize phi[eq=0] to 1 everywhere and 0 for other equations
	memset(phi, 0, (NEQNS * NCELLS) * sizeof(double));
	for (int i = 0; i < NCELLS; i++)
	{
		phi[IDX(i, 0, NCELLS)] = 1;
	}




	// ------ Write output file --------
	err = write_vtk_output("output_file.vtk", &nodes, &cells, &NPOINTS, &NCELLS,
		&CELL_LIST_SIZE, &phi);


	// Release Allocated Memory for grid
	free(nodes);
	free(cells);

	// Release conservative scalars memory
	free(phi);

	

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

/*---------------------------------------------------------------------------
* Write data function
----------------------------------------------------------------------------*/
int write_vtk_output(const char* out_filename,	node** nodes,	cell** cells,
	int* NPOINTS,	int* NCELLS,	int* CELL_LIST_SIZE,	double** phi)
{
	//Open file for writing
	FILE* fp = fopen(out_filename, "w");

	if (!fp)
	{
		perror("Error writing data file\n");
		return 1;
	}

	// Header section
	fprintf(fp, "# vtk DataFile Version 3.0\n");
	fprintf(fp, "Jack's first unstructured vtk datafile\n");
	fprintf(fp, "ASCII\n");

	// Dataset type definition
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

	// Write points
	fprintf(fp, "POINTS %d double\n", *NPOINTS);
	
	for (int i = 0; i < *NPOINTS; i++) 
	{
		fprintf(fp, "%.15f %.15f %.15f\n", (*nodes)[i].x, (*nodes)[i].y, (*nodes)[i].z);
	}

	// Write Cells
	fprintf(fp, "CELLS %d %d\n", *NCELLS, *CELL_LIST_SIZE);
	
	for (int i = 0; i < *NCELLS; i++)
	{
		// Get number of nodes for current cell
		int num_nodes = (*cells)[i].num_nodes;
		fprintf(fp, "%d ", num_nodes);

		for (int point = 0; point < num_nodes; point++)
		{
			if (point == (num_nodes - 1))
			{
				fprintf(fp, "%d\n", (*cells)[i].node_ids[point]);
			}
			else
			{
				fprintf(fp, "%d ", (*cells)[i].node_ids[point]);
			}
			
		}
	}

	// Write Cell Types
	fprintf(fp, "CELL_TYPES %d \n", *NCELLS);

	for (int i = 0; i < *NCELLS; i++)
	{
		fprintf(fp, "%d\n", (*cells)[i].type);
	}

	// Write Scalar data 
	fprintf(fp, "CELL_DATA %d\n", *NCELLS);

	fprintf(fp, "SCALARS phi[%d] double 1\n",0);
	fprintf(fp, "LOOKUP_TABLE default\n");

	for (int i = 0; i < *NCELLS; i++)
	{
		fprintf(fp, "%.15f\n", (*phi)[IDX(i, 0, *NCELLS)]);
	}

	fclose(fp);

	return 0;
}