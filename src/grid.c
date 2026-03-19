/*-------------------------------------------------------
* grid.c - Implementation of grid data structures and functions
* -------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "grid.h"

int get_num_faces(int vtk_type)
{
	switch (vtk_type)
	{
		//Degenerate Cells 
	case VTK_EMPTY_CELL:
	case VTK_VERTEX:
	case VTK_POLY_VERTEX:
	case VTK_LINE:
	case VTK_POLY_LINE:
		return 0;

		// 2D Cells 
	case VTK_TRIANGLE:
		return 3;

	case VTK_PIXEL:
		return 4;

	case VTK_QUAD:
		return 4;

		// 3D Cells

	case VTK_TETRA:
		return 4;
	case VTK_VOXEL:
		return 6;
	case VTK_HEXAHEDRON:
		return 6;
	case VTK_WEDGE:
		return 5;
	case VTK_PYRAMID:
		return 5;
	case VTK_PENTAGONAL_PRISM:
		return 7;
	case VTK_HEXAGONAL_PRISM:
		return 8;

	default:
		fprintf(stderr, "Unsupported VTK cell type %d\n", vtk_type);
		exit(1);
	}
}

// Read grid .vtk file
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
			for (int i = 0; i < *NCELLS; i++)
			{
				cells[i].node_ids = NULL; // Initialize pointer to zero
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

			/*
			for (int i = 0; i < 8; i++)
			{
				cells[cidx].node_ids[i] = -1; // Initialize node IDs to -1
			}*/

			if (sscanf(line, "%d %d %d %d %d %d %d %d %d", &num_nodes, &node_ids[0], &node_ids[1], &node_ids[2], &node_ids[3], &node_ids[4], &node_ids[5], &node_ids[6], &node_ids[7]) >= 2)
			{

				// Store node IDs for the cell
				cells[cidx].num_nodes = num_nodes; // Store number of nodes in the cell

				// allocate node_ids 
				cells[cidx].node_ids = malloc(num_nodes * sizeof(int));

				for (int i = 0; i < num_nodes; i++)
				{
					cells[cidx].node_ids[i] = node_ids[i];

				}

				// Store cell ID
				cells[cidx].id = cidx;							
								 
				// Calculate centroid, Cell Volume, Fill out face information also
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
	if (nodes && pidx != *NPOINTS) { free_grid(nodes, cells, cidx); return 11; }
	if (cells && cidx != *NCELLS) { free_grid(nodes, cells, cidx); return 12; }

	// Set output pointers
	*nodes_out = nodes;
	*cells_out = cells;

	return 0; // Success

}

void free_grid(node* nodes, cell* cells, int ncells_allocated)
{
	free(nodes);

	if (cells)
	{
		for (int i = 0; i < ncells_allocated; ++i)
		{
			free(cells[i].node_ids);
		}
		free(cells);
	}
}

// Function to calculate cell centroid, volume, face information, and other geometric properties
int calculate_cell_centroid_and_volume(node* nodes, cell* cells, int *NCELLS)
{
	// Loop over cells
	//	calculate geometric center of cell using 6.21 
	//	

	for (int i = 0; i < *NCELLS; i++)
	{
		cell* c = &cells[i];
		
		// Calculate number of faces for given cell
		c->num_faces = get_num_faces(c->type);

		// Calculate geometric center
		double x_G[3] = { 0.0, 0.0, 0.0 };

		for (int j = 0; j < c->num_nodes; j++)
		{
			int node_id = c->node_ids[j];
			x_G[0] += nodes[node_id].x;
			x_G[1] += nodes[node_id].y;
			x_G[2] += nodes[node_id].z;
		}

		x_G[0] /= c->num_nodes;
		x_G[1] /= c->num_nodes;
		x_G[2] /= c->num_nodes;
		
		cblas_dscal(3, 1.0 / c->num_nodes, x_G, 1);
				
		// Get number of sub angles
		int num_sub_angles = (c->num_faces)*(c->num_faces - 1)*(c->num_faces - 2)/6;

		// Calculate total cell volume by summing sub-volumes of tetrahedra formed by cell centroid and each face

		


	}

	return 0;
}