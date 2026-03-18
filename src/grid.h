/*-------------------------------------------------------
* grid.h - Header file for grid data structures and functions
* -------------------------------------------------------*/

#ifndef GRID_H
#define GRID_H

typedef struct node {
	double x, y, z; // Node coordinates
	int id; // Node ID
} node;


// Cell Types
enum {
	// Linear cells
	VTK_EMPTY_CELL = 0,
	VTK_VERTEX = 1,
	VTK_POLY_VERTEX = 2,
	VTK_LINE = 3,
	VTK_POLY_LINE = 4,
	VTK_TRIANGLE = 5,
	VTK_TRIANGLE_STRIP = 6,
	VTK_POLYGON = 7,
	VTK_PIXEL = 8,
	VTK_QUAD = 9,
	VTK_TETRA = 10,
	VTK_VOXEL = 11,
	VTK_HEXAHEDRON = 12,
	VTK_WEDGE = 13,
	VTK_PYRAMID = 14,
	VTK_PENTAGONAL_PRISM = 15,
	VTK_HEXAGONAL_PRISM = 16,
};

typedef struct cell {
	int* node_ids; // IDs of the 8 nodes that form the cell
	int type; // Cell type (e.g., hexahedron, tetrahedron)
	int id; // Cell ID
	int num_nodes; // Number of nodes in the cell (e.g., 8 for hexahedron, 4 for tetrahedron)

	double xc, yc, zc; // Cell centroid
	double volume;
	int* face_ids;
	int num_faces;
} cell;

typedef struct face {
	int* node_ids;
	int num_nodes;
	int owner;
	int neighbor;

	double xc, yc, zc; //face centroid
	double sx, sy, sz; //Face area vector
}face;

int read_grid(const char* filename,
	node** nodes_out,
	cell** cells_out,
	int* NPOINTS,
	int* NCELLS,
	int* CELL_LIST_SIZE);

void free_grid(node* nodes,
	cell* cells,
	int ncells_allocated);

int get_num_faces(int vtk_type);

#endif