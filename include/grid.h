/*-------------------------------------------------------
* grid.h - Header file for grid data structures and functions
* -------------------------------------------------------*/

#ifndef GRID_H
#define GRID_H

#include <stdbool.h>
#include "mkl.h"	
#include "math_helpers.h"


// Boundary Types
typedef enum boundaryType{
	Dirichlet = 0,
	Neumann = 1,
	Robin = 2,
}boundaryType;

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

typedef struct boundary {
	int id; // Surface ID
	int endpoints[2]; // Endpoints of the surface (if applicable)
	int num_faces; // Number of faces in the surface
	int* face_ids; // Array of face IDs that define the surface
	boundaryType type; // Boundary condition type (Dirichlet, Neumann, Robin)
} boundary;


typedef struct node {
	double x, y, z; // Node coordinates
	int id; // Node ID
} node;

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
	int id; 
	int num_nodes;
	int owner; //owner cell ID
	int neighbor; //neighbor cell ID (uses cell id of degenerate face/edge cell for boundary faces)
	bool boundary_face; //flag to indicate if this is a boundary face

	double xc, yc, zc; //face centroid
	double Sx, Sy, Sz; //Face area vector
	double Ex, Ey, Ez; // orthogonal contribution
	double Tx, Ty, Tz; // tangential contribution
}face;

int read_grid(const char* filename, //grid filename INPUT
	node** nodes_out, 
	cell** cells_out,
	int* NPOINTS,
	int* NCELLS,
	int* CELL_LIST_SIZE,
	int* MAX_FACES,
	int* NDEGEN_CELLS);

void free_grid(node* nodes,
	cell* cells,
	face* faces,
	int ncells_allocated,
	int nfaces_allocated);

int get_num_faces(int vtk_type);

int calculate_cell_centroid_and_vol(cell* c, node* nodes);

int build_interior_face(cell* c, face* faces, node* nodes,cell* cells, int k, int* fidx);

int build_boundary_face(cell* c, face* faces, node* nodes, cell* cells, int k, int* fidx);

int build_face(cell* c, face* faces, node* nodes, cell* cells, int k, int* fidx);

int build_boundary(boundary* b, int id, int* endpoints, boundaryType type, node* nodes, face* faces, int* NFACES);

int build_faces_and_cells(node* nodes,
	cell* cells,
	int* NCELLS,
	int* MAX_FACES,
	int* NFACES, // Number of faces in grid OUTPUT
	face** faces_out); // Face connectivity OUTPUT

// Calculate Face Centroid and Area Vector 
int calculate_FC_AV(face* f, cell* c_owner, cell* c_neighbor, node* nodes, int* node_ids);

int comp(const void* a, const void* b); //Comparator function for ascending order sorting qsort 



#endif