#ifndef CFD_H
#define CFD_H

typedef struct node {
	double x, y, z; // Node coordinates
	int id; // Node ID
} node;

typedef struct cell {
	int node_ids[8]; // IDs of the 8 nodes that form the cell
	int type; // Cell type (e.g., hexahedron, tetrahedron)
	int id; // Cell ID
	int num_nodes; // Number of nodes in the cell (e.g., 8 for hexahedron, 4 for tetrahedron)
} cell;

int read_grid(const char* filename,
	node** nodes_out,
	cell** cells_out,
	int* NPOINTS,
	int* NCELLS,
	int* CELL_LIST_SIZE);


#endif // !