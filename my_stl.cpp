#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include "admesh/stl.h"
#include <string.h>


int my_stl_compare_function(stl_hash_edge *edge_a, stl_hash_edge *edge_b) {
	if (edge_a->facet_number == edge_b->facet_number) {
		return 1;			/* Don't match edges of the same facet */
	}
	else {
		return memcmp(edge_a, edge_b, SIZEOF_EDGE_SORT);
		/*if (!memcmp(&edge_a->key[0], &edge_b->key[0], 3 * sizeof(unsigned)) && !memcmp(&edge_a->key[1], &edge_b->key[1], 3 * sizeof(unsigned))) {
		return 0;
		}
		if (!memcmp(&edge_a->key[0], &edge_b->key[1], 3 * sizeof(unsigned)) && !memcmp(&edge_a->key[1], &edge_b->key[0], 3 * sizeof(unsigned))) {
		printf("Reverse\n");
		return 0;
		}
		return 1;*/
	}
}

int my_stl_get_hash_for_edge(int M, stl_hash_edge *edge) {
	return ((edge->key[0] / 23 + edge->key[1] / 19 + edge->key[2] / 17
		+ edge->key[3] / 13 + edge->key[4] / 11 + edge->key[5] / 7) % M);
}

void my_insert_hash_edge(stl_file *stl, stl_hash_edge edge) {
	stl_hash_edge *link;
	stl_hash_edge *new_edge;
	int chain_number;

	chain_number = my_stl_get_hash_for_edge(stl->M, &edge);

	link = stl->heads[chain_number];

	if (link == stl->tail) {
		/* This list doesn't have any edges currently in it.  Add this one. */
		new_edge = (stl_hash_edge*)malloc(sizeof(stl_hash_edge));
		if (new_edge == NULL) perror("my_insert_hash_edge");
		stl->stats.malloced++;
		*new_edge = edge;
		new_edge->next = stl->tail;
		stl->heads[chain_number] = new_edge;
		return;
	}
	else  if (!my_stl_compare_function(&edge, link)) {
		// simply return, since we want a list of edges
		stl->stats.collisions++;
		return;
	}
	else {
		/* Continue through the rest of the list */
		for (;;) {
			if (link->next == stl->tail) {
				/* This is the last item in the list. Insert a new edge. */
				new_edge = (stl_hash_edge*)malloc(sizeof(stl_hash_edge));
				if (new_edge == NULL) perror("insert_hash_edge");
				stl->stats.malloced++;
				*new_edge = edge;
				new_edge->next = stl->tail;
				link->next = new_edge;
				return;
			}
			else  if (!my_stl_compare_function(&edge, link->next)) {
				// simply return, since we want a list of edges
				stl->stats.collisions++;
				return;
			}
			else {
				/* This is not a match.  Go to the next link */
				link = link->next;
			}
		}
	}
}


void my_stl_initialize_facet_check_exact(stl_file *stl) {
	int i;

	if (stl->error) return;

	stl->stats.malloced = 0;
	stl->stats.freed = 0;
	stl->stats.collisions = 0;


	stl->M = 81397;

	stl->heads = (stl_hash_edge**)calloc(stl->M, sizeof(*stl->heads));
	if (stl->heads == NULL) perror("stl_initialize_facet_check_exact");

	stl->tail = (stl_hash_edge*)malloc(sizeof(stl_hash_edge));
	if (stl->tail == NULL) perror("stl_initialize_facet_check_exact");

	stl->tail->next = stl->tail;

	for (i = 0; i < stl->M; i++) {
		stl->heads[i] = stl->tail;
	}
}

void my_stl_load_edge_exact(stl_file *stl, stl_hash_edge *edge,
	stl_vertex *a, stl_vertex *b) {

	float diff_x;
	float diff_y;
	float diff_z;
	float max_diff;

	if (stl->error) return;

	diff_x = ABS(a->x - b->x);
	diff_y = ABS(a->y - b->y);
	diff_z = ABS(a->z - b->z);
	max_diff = STL_MAX(diff_x, diff_y);
	max_diff = STL_MAX(diff_z, max_diff);
	stl->stats.shortest_edge = STL_MIN(max_diff, stl->stats.shortest_edge);

	if (diff_x == max_diff) {
		if (a->x > b->x) {
			memcpy(&edge->key[0], a, sizeof(stl_vertex));
			memcpy(&edge->key[3], b, sizeof(stl_vertex));
		}
		else {
			memcpy(&edge->key[0], b, sizeof(stl_vertex));
			memcpy(&edge->key[3], a, sizeof(stl_vertex));
			edge->which_edge += 3; /* this edge is loaded backwards */
		}
	}
	else if (diff_y == max_diff) {
		if (a->y > b->y) {
			memcpy(&edge->key[0], a, sizeof(stl_vertex));
			memcpy(&edge->key[3], b, sizeof(stl_vertex));
		}
		else {
			memcpy(&edge->key[0], b, sizeof(stl_vertex));
			memcpy(&edge->key[3], a, sizeof(stl_vertex));
			edge->which_edge += 3; /* this edge is loaded backwards */
		}
	}
	else {
		if (a->z > b->z) {
			memcpy(&edge->key[0], a, sizeof(stl_vertex));
			memcpy(&edge->key[3], b, sizeof(stl_vertex));
		}
		else {
			memcpy(&edge->key[0], b, sizeof(stl_vertex));
			memcpy(&edge->key[3], a, sizeof(stl_vertex));
			edge->which_edge += 3; /* this edge is loaded backwards */
		}
	}
}

void my_stl_free_edges(stl_file *stl) {
	int i;
	stl_hash_edge *temp;

	if (stl->error) return;

	if (stl->stats.malloced != stl->stats.freed) {
		for (i = 0; i < stl->M; i++) {
			for (temp = stl->heads[i]; stl->heads[i] != stl->tail;
				temp = stl->heads[i]) {
				stl->heads[i] = stl->heads[i]->next;
				free(temp);
				stl->stats.freed++;
			}
		}
	}
	free(stl->heads);
	free(stl->tail);
}