/*-------------------------------------------------------
* grid.c - Implementation of grid data structures and functions
* -------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "math_helpers.h"
#include "grid.h"

int get_num_faces(int vtk_type)
{
	switch (vtk_type)
	{
		//Degenerate Cells 
	case VTK_EMPTY_CELL:
	case VTK_VERTEX:
	case VTK_POLY_VERTEX:
		return 0;
	case VTK_LINE:
	case VTK_POLY_LINE:
		return 1;

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
int read_grid(const char* filename, node** nodes_out, cell** cells_out, int* NPOINTS, int* NCELLS, int* CELL_LIST_SIZE, int* MAX_FACES, int* NDEGEN_CELLS)
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
	bool cell_entity_section = false;

	int pidx = 0; // Node index
	int cidx = 0; // Cell index
	int ctidx = 0; // Cell type index
	int ceidx = 0; // Cell entity index

	// Create Null Pointers for nodes and cells, will allocate memory after reading number of points and cells
	node* nodes = NULL;
	cell* cells = NULL;

	// Initialize MAX_FACES
	*MAX_FACES = 0;

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
			cell_entity_section = false;

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
			cell_entity_section = false;

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
				// Initialize pointers to zero
				cells[i].node_ids = NULL;
				cells[i].face_ids = NULL;
			}

			cidx = 0; // Reset cell index for reading cell data
			continue;
		}

		if (sscanf(line, "CELL_TYPES %d", NCELLS) == 1)
		{
			cell_types_section = true;
			points_section = false;
			cells_section = false;
			cell_entity_section = false;

			ctidx = 0; // Reset cell type index for reading cell type data
			continue;
		}

		if (sscanf(line, "CELL_DATA %d", NCELLS) == 1)
		{
			cell_entity_section = true;
			cell_types_section = false;
			points_section = false;
			cells_section = false;

			ceidx = 0; // Reset cell entity index for reading cell entity data
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

			if (cell_type < 5) // If Cell is a no volume cell, increment degenerate cell count
			{
				(*NDEGEN_CELLS)++;
			}

			*MAX_FACES += get_num_faces(cell_type);

			ctidx++;
		}
		else if (cell_entity_section)
		{
			int entity_id;
			if (sscanf(line, "%d", &entity_id) == 1)
			{
				cells[ceidx].entity_id = entity_id;
			}
			else if (strcmp(line, "SCALARS CellEntityIds int 1\n") == 0);
			else if (strcmp(line, "LOOKUP_TABLE default\n") == 0);
			else
			{
				fprintf(stderr, "Error reading CellEntityIds data\n");
				return 1;
			}

			ceidx++;
		}
	}

	// close file
	fclose(fp);

	//Release memory on error 
	if (nodes && pidx != *NPOINTS) { free_grid(nodes, cells, NULL, cidx, 0); return 11; }
	if (cells && cidx != *NCELLS) { free_grid(nodes, cells, NULL, cidx, 0); return 12; }

	// Set output pointers
	*nodes_out = nodes;
	*cells_out = cells;

	return 0; // Success

}

// Free Grid related variables from memory function 
void free_grid(node* nodes, cell* cells, face* faces, int ncells_allocated, int nfaces_allocated)
{
	free(nodes);

	if (cells)
	{
		for (int i = 0; i < ncells_allocated; ++i)
		{
			free(cells[i].node_ids);
			free(cells[i].face_ids);
		}
		free(cells);
	}

	if (faces)
	{
		for (int i = 0; i < nfaces_allocated; ++i)
		{
			free(faces[i].node_ids);
		}
		free(faces);
	}
}

// Function to calculate cell centroid, volume, face information, and other geometric properties
int build_faces_and_cells(node* nodes, cell* cells, int* NCELLS, int* MAX_FACES, int* NFACES, face** faces_out)
{
	// Loop over cells
	//	calculate geometric center of cell using 6.21 
	//	

	// Allocate faces based on worse case scenario
	face* faces = malloc((size_t)(*MAX_FACES) * sizeof(face));

	if (faces == NULL)
	{
		fprintf(stderr, "Error: memory allocation failed for faces array \n");
		return 2;
	}

	// Initialize node_ids pointer to zero for face 
	for (int i = 0; i < *MAX_FACES; i++)
	{
		faces[i].node_ids = NULL;
	}

	*NFACES = 0; // Initialize number of faces
	//int mfidx = 0; // Face array index (indexs the max size array)
	int fidx = 0; // actual face index (actual face id)

	// Loop over all cells
	for (int i = 0; i < *NCELLS; i++)
	{
		cell* c = &cells[i];
		
		// Check for degenerate cell
		/*
		if (c->type < 5) // If Cell is a no volume cell
		{
			c->volume = 0.0;
			continue;
		}*/

		// Calculate number of faces for given cell and allocate face ids
		c->num_faces = get_num_faces(c->type);
		c->face_ids = malloc(c->num_faces * sizeof(int));
		if (c->face_ids == NULL)
		{
			fprintf(stderr, "Error allocating face id memory \n");
			// clean up anything already allocated
			for (int k = 0; k < i; ++k)
			{
				free(cells[k].face_ids);
			}
			return 2;
		}

		// Calculate total cell volume by summing sub-volumes triangles formed by cell centroid and each face
		// Also allocate face information

		/*----------------- Compute Cell Volume and Centroid ----------------------*/
		int err = calculate_cell_centroid_and_vol(c, nodes);

		for (int k = 0; k < c->num_faces; k++)
		{
			/*-------------Compute Face Connectivity---------------------*/
			//err = build_interior_face(c, faces, nodes, k, &fidx);
			err = build_face(c, faces, nodes, cells, k, &fidx);
			if (err != 0)
			{
				free(faces);
				fprintf(stderr, "Error building faces \n");
				return 2;

			}	

		}	
		
	}
	
	// Set the number of faces and face array 
	*NFACES = fidx;

	// Reallocate faces and assign it to faces_out
	face* tmp = realloc(faces, (size_t)fidx * sizeof(face));
	if (tmp != NULL)
	{
		faces = tmp;
	}

	*faces_out = faces; 
	return 0;
}

int calculate_cell_centroid_and_vol(cell* c, node* nodes)
{
	// Calculate geometric center
	double x_G[3] = { 0.0, 0.0, 0.0 };

	for (int j = 0; j < c->num_nodes; j++)
	{
		int node_id = c->node_ids[j];
		x_G[0] += nodes[node_id].x;
		x_G[1] += nodes[node_id].y;
		x_G[2] += nodes[node_id].z;
	}

	cblas_dscal(3, 1.0 / c->num_nodes, x_G, 1);

	if (c->type < 5) // If Cell is a degenerate cell, centroid is geometric center and volume is zero
	{
		c->volume = 0.0;
		c->xc = x_G[0];
		c->yc = x_G[1];
		c->zc = x_G[2];

		return 0; //Function ends here for degenerate cells
	} 

	// Initialize cell volume and cell centroid coordinates to zero
		c->volume = 0.0;
		c->xc = 0.0;
		c->yc = 0.0;
		c->zc = 0.0;

	for (int k = 0; k < c->num_faces; k++)
	{
		/*----------------- Compute Cell Volume and Centroid ----------------------*/
		// Get node IDs for the current sub-triangle (Also node ids of the first face)
		int node_id1 = c->node_ids[k % c->num_faces];
		int node_id2 = c->node_ids[(k + 1) % c->num_faces];

		// Get coordinates of the nodes
		double r1[3] = { nodes[node_id1].x, nodes[node_id1].y, nodes[node_id1].z };
		double r2[3] = { nodes[node_id2].x, nodes[node_id2].y, nodes[node_id2].z };
		double r3[3] = { x_G[0], x_G[1], x_G[2] };

		// Calculate geometric center/centroid of subtriangle
		double x_CE_t = (r1[0] + r2[0] + r3[0]) / 3.0;
		double y_CE_t = (r1[1] + r2[1] + r3[1]) / 3.0;
		double z_CE_t = (r1[2] + r2[2] + r3[2]) / 3.0;


		// Compute area of subtriangle formed by r1, r2, and r3 using cross product
		double v1[3];
		double v2[3];

		cblas_dcopy(3, r2, 1, v1, 1); // v1 = r2
		cblas_dcopy(3, r3, 1, v2, 1); // v2 = r3

		cblas_daxpy(3, -1.0, r1, 1, v1, 1); // v1 = r2 - r1
		cblas_daxpy(3, -1.0, r1, 1, v2, 1); // v2 = r3 - r1

		// Calculate volume of the triangle formed by r1, r2, and r3
		double St[3];
		cross_prod(v1, v2, St); // Calculate cross product of v1 and v2to get the area vector of the triangle

		cblas_dscal(3, 0.5, St, 1); // Scale the area vector by 0.5 to get the area of the triangle
		double St_mag;
		magnitude(St, &St_mag); // Calculate the magnitude of the area vector to get the area of the triangle

		// Add Triangle volume/area to total cell volume
		c->volume += St_mag;

		// Add sub triangle contrubtion to cell centroid calculation
		c->xc += x_CE_t * St_mag;
		c->yc += y_CE_t * St_mag;
		c->zc += z_CE_t * St_mag;
	}

	// Divide by total cell volume to get final cell centroid
	c->xc /= c->volume;
	c->yc /= c->volume;
	c->zc /= c->volume;

	return 0;
}

int build_face(cell* c, face* faces, node* nodes, cell* cells, int k, int* fidx)
{
	faces[*fidx].boundary_face = false; // Initialize boundary face flag to false, will be set to true for boundary faces in build_boundary_face
	faces[*fidx].boundary_id = -1; //Initialize boundary id to -1, will be set later when boundaries are built 

	if (c->type < 5) //Check for degenerate cell
	{
		return build_boundary_face(c, faces, nodes, cells, k, fidx); // Build boundary face for degenerate cell
	}
	else
	{
		return build_interior_face(c, faces, nodes, cells,  k, fidx); // Build interior face for normal cell
	}

}

int build_boundary_face(cell* c, face* faces, node* nodes, cell* cells, int k, int* fidx)
{
	if (c->type < 2) // Check for single node degenerate cell (corner)
	{
		return 0; // Do not build face for vertex like degenerate cells
	}

	// Get number of face nodes (always 2 for 2d problem)
	int num_nodes = 2;

	// Otherwise, its a line degenerate cell. Assign the cell id as the neighbor cell of the face
	faces[*fidx].neighbor = c->id;
	int node_id1 = c->node_ids[0];
	int node_id2 = c->node_ids[1];
	int node_ids[2] = { node_id1, node_id2 };
	qsort(node_ids, 2, sizeof(int), comp); //Sort node ids in ascending order

	faces[*fidx].num_nodes = num_nodes;
	faces[*fidx].node_ids = malloc(num_nodes * sizeof(int));
	if (faces[*fidx].node_ids == NULL)
	{
		fprintf(stderr, "Error allocating node ararys for faces\n");

		for (int alo_fidx = 0; alo_fidx < *fidx; alo_fidx++)
		{
			free(faces[alo_fidx].node_ids);
		}
		//free(faces); Should free faces in caller

		return 2;
	}

	// Set face data 
	faces[*fidx].node_ids[0] = node_ids[0];
	faces[*fidx].node_ids[1] = node_ids[1];

	faces[*fidx].boundary_face = true;
	faces[*fidx].id = *fidx;
	c->face_ids[k] = *fidx; // add face to cell

	(*fidx)++; // Increment face index for next face

	// All other properties will be assigned in build_interior_face when the face is first created, since the face will only be seen once for boundary faces and will be created as a new face with the current cell as the owner. This will only work for degenerate cells with > 2 nodes. Will not work for vertex like cells
	// Vertex cells get -1 as the neighbor

	return 0;
}
// Make a build face function that checks if the cell is degenerate or not, then 
// do build_boundary face if degenerate which computes the surface vector and 
// assigns the degenrate cell as the neighbor cell. This will only work for 
// degenerate cells with > 2 nodes. Will not work for vertex like cells 
// if its a normal cell, use build interior face


int build_interior_face(cell* c, face* faces, node* nodes, cell* cells, int k, int* fidx)
{
	// Get number of face nodes (always 2 for 2d problem)
	int num_nodes = 2;

	// Calculate number of new faces 
	// Get node IDs for the current sub-triangle (Also node ids of the first face)
	int node_id1 = c->node_ids[k % c->num_faces];
	int node_id2 = c->node_ids[(k + 1) % c->num_faces];
	int node_ids[2] = { node_id1, node_id2 };
	qsort(node_ids, 2, sizeof(int), comp); //Sort node ids in ascending order

	// Flags for if face is found
	bool oldFaceFlag = false;
	int oldFaceidx = 0;

	// Loop over existing faces to check if this one has been found, is this face already in the list of faces?
	for (int l = 0; l < *fidx; l++)
	{
		if ((faces[l].node_ids[0] == node_ids[0])
			&& (faces[l].node_ids[1] == node_ids[1]))
		{
			// Does this face already have an neighbor cell?
			// Are the nodes of this face the same as the nodes of the current face (regardless of order)? If so, then this is not a new face, and the current cell is the neighbor cell of this face. This will not trigger for boundary faces since they will only be seen once and will have neighbor -1.
			// current cell is neighbor cell
			// leave owner cell (should already be defined

			oldFaceFlag = true;
			oldFaceidx = l;
			break;
		}
	}

	// Decide action based on if face is new
	if (oldFaceFlag)
	{
		// Old Face, current cell is neigbor cell (This will not trigger for boundary faces)
		// Everything else should be allocated on the first pass
		
		if (faces[oldFaceidx].boundary_face) // If it is a boundary face, it will have been seen as a degenerate cell. Assign the owner cell as the non-degnerate cell and compute its face area vector
		{
			faces[oldFaceidx].owner = c->id;

			cell* c_neighbor = &cells[faces[oldFaceidx].neighbor]; // Get the neighbor cell (which is the degenerate cell that was used to build the boundary face)
			calculate_FC_AV(&faces[oldFaceidx], c, c_neighbor, nodes, node_ids);
		}
		else
		{
			faces[oldFaceidx].neighbor = c->id;

			cell* c_owner = &cells[faces[oldFaceidx].owner]; // Get the owner cell (which is the first cell that was used to build the face)
			calculate_FC_AV(&faces[oldFaceidx], c_owner, c, nodes, node_ids);
		}

		// Add the face id to the cell. old face index is added to the cell
		c->face_ids[k] = faces[oldFaceidx].id;

		// Reset Flags 
		oldFaceFlag = false;
		oldFaceidx = 0;
	}
	else
	{
		// New Face, current cell is owner
		// No neighbor yet
		// ID is the newest face 
		// Boundary faces will not be seen twice and will have neighbor -1.
		// 
		// Allocate node_ids(should only be done for a new face)
		faces[*fidx].num_nodes = num_nodes;
		faces[*fidx].node_ids = malloc(num_nodes * sizeof(int));
		if (faces[*fidx].node_ids == NULL)
		{
			fprintf(stderr, "Error allocating node ararys for faces\n");

			for (int alo_fidx = 0; alo_fidx < *fidx; alo_fidx++)
			{
				free(faces[alo_fidx].node_ids);
			}
			//free(faces); Should free faces in caller

			return 2;
		}

		// Set face data 
		faces[*fidx].node_ids[0] = node_ids[0];
		faces[*fidx].node_ids[1] = node_ids[1];

		faces[*fidx].owner = c->id;

		// Only make neighbor -1 if the neighbor id has not been already set. 
		// A face on a degenerate cell will already be seen in build_boundary_face, 
		// but the neighbor will be assigned to the degenerate cell id, so we do not want to overwrite that with -1.
		faces[*fidx].neighbor = -1;
		//if (!faces[*fidx].neighbor)
		//{
		//	faces[*fidx].neighbor = -1;
		//}
		
		faces[*fidx].id = *fidx;

		//// Add the face id to the cell. New face index is added to the cell
		c->face_ids[k] = *fidx;

		//increment face counter
		(*fidx)++;

	}
	return 0;
}

int calculate_FC_AV(face* f, cell* c_owner, cell* c_neighbor, node* nodes, int* node_ids)
{
	// Face centroid
	f->xc = (nodes[node_ids[0]].x + nodes[node_ids[1]].x) / 2;
	f->yc = (nodes[node_ids[0]].y + nodes[node_ids[1]].y) / 2;
	f->zc = (nodes[node_ids[0]].z + nodes[node_ids[1]].z) / 2;

	// Face Surface Vector components
	double dx = nodes[node_ids[1]].x - nodes[node_ids[0]].x;
	double dy = nodes[node_ids[1]].y - nodes[node_ids[0]].y;

	// Create tangent surface vector
	double T[3] = { dx, dy };

	// Rotate the vector 90 degrees to get the normal vector candidate
	double Sf[3] = { dy, -dx, 0 }; // Surface area of face is the cross product of v1 and v3, 1/2((r2 - r1) x (0 - r1))

	// Check if the vector is point in the right direction (out from owner cell)
	if (((f->xc - c_owner->xc) * Sf[0] + (f->yc - c_owner->yc) * Sf[1]) < 0)
	{
		// If the dot product is negative, the face normal is pointing inward, so we need to flip it
		Sf[0] = -Sf[0];
		Sf[1] = -Sf[1];
	}

	// Assign face vector to faces
	f->Sx = Sf[0];
	f->Sy = Sf[1];
	f->Sz = Sf[2];

	double Sf_mag = 0;
	magnitude(Sf, &Sf_mag);

	//Use Orthogonal Correction Approach
	double n_hat[3] = { f->Sx/Sf_mag, f->Sy/Sf_mag, f->Sz/Sf_mag }; // face unit normal vector

	// Build e vector
	// Calculate vector from owner cell centroid to neighbor cell centroid
	double r_ON[3] = { c_neighbor->xc - c_owner->xc, c_neighbor->yc - c_owner->yc, c_neighbor->zc - c_owner->zc };

	double r_ON_mag;
	magnitude(r_ON, &r_ON_mag);
	double e[3] = { r_ON[0] / r_ON_mag, r_ON[1] / r_ON_mag, r_ON[2] / r_ON_mag };

	// Orthogonal contribution vector (use orthogonal correction approach)
	f->Ex = Sf_mag*e[0];
	f->Ey = Sf_mag*e[1];
	f->Ez = Sf_mag*e[2];

	// Cross Diffusion contribution vector (use normal correction) eq 8.71
	f->Tx = Sf[0] - f->Ex;
	f->Ty = Sf[1] - f->Ey;
	f->Tz = Sf[2] - f->Ez;

	return 0;
}

int build_boundary(boundary* b, int id, int* endpoints, boundaryType type, boundaryData bData, node* nodes, face* faces, int* NFACES)
{

	int initial_capacity = 10; // Initial capacity for face_ids array

	b->id = id;
	b->endpoints[0] = endpoints[0];
	b->endpoints[1] = endpoints[1];
	b->type = type;

	// initialize face ids
	b->face_ids = NULL;
	b->face_ids = malloc(initial_capacity * sizeof(int)); // Will realloc later based on number of faces found
	b->num_faces = 0; 

	// Boundary line
	double x0 = nodes[endpoints[0]].x;
	double x1 = nodes[endpoints[1]].x;

	double y0 = nodes[endpoints[0]].y;
	double y1 = nodes[endpoints[1]].y;


	if (x0 == x1) // Vertical Line
	{
		// Loop over faces and find faces that have x coordinate between x0 and x1 and y coordinate between y0 and y1
		for (int i = 0; i < *NFACES; i++)
		{
			if (faces[i].xc >= x0 - 1e-6 && faces[i].xc <= x0 + 1e-6)
			{
				// Check if we need to realloc face_ids array
				if (b->num_faces >= initial_capacity)
				{
					initial_capacity *= 2;
					int* tmp = realloc(b->face_ids, initial_capacity * sizeof(int));
					if (tmp == NULL)
					{
						fprintf(stderr, "Error reallocating memory for boundary face ids\n");
						free(b->face_ids);
						return 2;
					}
					b->face_ids = tmp;
				}
				b->face_ids[b->num_faces] = faces[i].id;
				faces[i].boundary_id = b->id; // Set the boundary id for the face
				b->num_faces++;
			}
		}

	}
	else if (y0 == y1) // Horizontal Line
	{
		// Loop over faces and find faces that have x coordinate between x0 and x1 and y coordinate between y0 and y1
		for (int i = 0; i < *NFACES; i++)
		{
			if (faces[i].yc >= y0 - 1e-6 && faces[i].yc <= y0 + 1e-6)
			{
				// Check if we need to realloc face_ids array
				if (b->num_faces >= initial_capacity)
				{
					initial_capacity *= 2;
					int* tmp = realloc(b->face_ids, initial_capacity * sizeof(int));
					if (tmp == NULL)
					{
						fprintf(stderr, "Error reallocating memory for boundary face ids\n");
						free(b->face_ids);
						return 2;
					}
					b->face_ids = tmp;
				}
				b->face_ids[b->num_faces] = faces[i].id;
				faces[i].boundary_id = b->id; // Set the boundary id for the face
				b->num_faces++;
			}
		}
	}
	else
	{
		fprintf(stderr, "Error: Boundary line is not vertical or horizontal\n");
		return 1;
	}

	// Realloc face_ids array to the correct size based on the number of faces found
	b->face_ids = realloc(b->face_ids, b->num_faces * sizeof(int));
	if ( b->face_ids == NULL)
	{
		fprintf(stderr, "Error reallocating memory for boundary face ids\n");
		free(b->face_ids);
		return 2;
	}

	//set boundary data
	b->data = bData;
	return 0;
}

// Math Helper functions (maybe move to separate heade
// Comparison function, must return negative if *a is less than *b
// 0 if they are equal
// and positive if *a is greater than *b
int comp(const void* a, const void* b) {
	return (*(int*)a - *(int*)b);
}


