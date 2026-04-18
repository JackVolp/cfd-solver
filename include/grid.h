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


typedef union boundaryData {
	double q_b; // Neumann boundary condition value (flux)
	struct {
		double h_inf; // Robin boundary condition coefficient
		double phi_inf; // Robin boundary condition ambient temperature
	} robin;
	double phi_b; // Dirichlet boundary condition value (temperature)
} boundaryData;

typedef struct boundary {
	int id; // Surface ID
	int endpoints[2]; // Endpoints of the surface (if applicable)
	int num_faces; // Number of faces in the surface
	int* face_ids; // Array of face IDs that define the surface
	boundaryType type; // Boundary condition type (Dirichlet, Neumann, Robin)
	// The conditional information if from ai, i havent read about how unions/ shared memory stuff works
	boundaryData data; // Union to hold boundary condition data (e.g., value for Dirichlet, flux for Neumann, coefficients for Robin)
} boundary;

typedef struct node {
	double x, y, z; // Node coordinates
	int id; // Node ID
} node;

typedef struct cellEntity {
	int id; // Cell Entity id
	int* cell_ids; // Array of cell ids that belong to the entity
	int num_cells; // Number of cells in the entity
} cellEntity;

typedef struct cell {
	int* node_ids; // IDs of the 8 nodes that form the cell
	int type; // Cell type (e.g., hexahedron, tetrahedron)
	int id; // Cell ID
	int num_nodes; // Number of nodes in the cell (e.g., 8 for hexahedron, 4 for tetrahedron)

	double xc, yc, zc; // Cell centroid
	double volume;
	int* face_ids;
	int num_faces;
	int entity_id;
} cell;

typedef struct face {
	int* node_ids;
	int id; 
	int num_nodes;
	int owner; //owner cell ID
	int neighbor; //neighbor cell ID (uses cell id of degenerate face/edge cell for boundary faces)
	bool boundary_face; //flag to indicate if this is a boundary face
	int boundary_id; // If boundary face, store boundary id for easy access to boundary condition information. Otherwise can be set to -1 or some other invalid value.
	double xc, yc, zc; //face centroid
	double Sx, Sy, Sz; //Face area vector
	double Ex, Ey, Ez; // orthogonal contribution
	double Tx, Ty, Tz; // tangential contribution
}face;

int read_grid(const char* filename, //grid filename INPUT
	node** nodes_out,
	cell** cells_out,
	cellEntity** cellEntities_out,
	int* NPOINTS,
	int* NCELLS,
	int* CELL_LIST_SIZE,
	int* MAX_FACES,
	int* NDEGEN_CELLS,
	int* NENTITIES);

void free_cell_entities(cellEntity* cellEntities, int nentities);

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

int build_boundary(boundary* b, int id, int* endpoints, boundaryType type, boundaryData bData, node* nodes, face* faces, int* NFACES);

int build_boundary_entity(boundary* b, int id, boundaryType type, boundaryData bData, cellEntity* entity, face* faces, int* NFACES);

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