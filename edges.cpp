#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include "admesh/stl.h"
#include "FIBR3Dapp.hpp"


void find_feature_edges(stl_file *stl, std::vector<feature_edge> &feature_edges)
{
	int i, j, neigh;
	float costheta;
	std::set<int> seen_facets;


	// find feature edges by inspecting neighbors
	for (i = 0; i < stl->stats.number_of_facets; i++) {
		if (seen_facets.find(i) != seen_facets.end())continue;
		// an unseen facet
		seen_facets.insert(i);
		if ((double)(stl->facet_start[i].vertex[2].x - stl->facet_start[i].vertex[1].x)*(double)(stl->facet_start[i].vertex[0].y - stl->facet_start[i].vertex[1].y)
			- (double)(stl->facet_start[i].vertex[2].y - stl->facet_start[i].vertex[1].y)*(double)(stl->facet_start[i].vertex[0].x - stl->facet_start[i].vertex[1].x) > 0) {
			// ccw oriented
			printf("CCW oriented\n");
			for (j = 0; j < 3; j++) {
				if (seen_facets.find(stl->neighbors_start[i].neighbor[j]) != seen_facets.end())continue;
				// an unseen facet
				neigh = stl->neighbors_start[i].neighbor[j];
				costheta = acos(mymax(mymin(stl->facet_start[i].normal.x*stl->facet_start[neigh].normal.x +
					stl->facet_start[i].normal.y*stl->facet_start[neigh].normal.y +
					stl->facet_start[i].normal.z*stl->facet_start[neigh].normal.z, 1.0f), -1.0f));
				if (costheta > DEG_TO_RAD(5.0) && costheta < DEG_TO_RAD(175.0)) { // cos(10^o)
																				  //printf("Cos= %.6f -- Theta=%f\n", (double)costheta, acos((double)costheta));
																				  // we have an angle bigger than 10 degrees, so facets common edge is a feature edge
																				  // find which edge is common

					feature_edge f_edge;
					f_edge.facet1_number = i;
					f_edge.facet2_number = neigh;
					f_edge.facet1_neighbor = j;
					f_edge.p[0] = stl->facet_start[i].vertex[j];
					f_edge.p[1] = stl->facet_start[i].vertex[(j + 1) % 3];

					// check for convexity
					v3 external =
						v3(stl->facet_start[i].normal.y*stl->facet_start[neigh].normal.z -
							stl->facet_start[i].normal.z*stl->facet_start[neigh].normal.y,
							stl->facet_start[i].normal.z*stl->facet_start[neigh].normal.x -
							stl->facet_start[i].normal.x*stl->facet_start[neigh].normal.z,
							stl->facet_start[i].normal.x*stl->facet_start[neigh].normal.y -
							stl->facet_start[i].normal.y*stl->facet_start[neigh].normal.x);
					v3 ab;
					if (f_edge.p[1].z < f_edge.p[0].z) {
						ab = v3(f_edge.p[0].x - f_edge.p[1].x, f_edge.p[0].y - f_edge.p[1].y, f_edge.p[0].z - f_edge.p[1].z);
					}
					else {
						ab = v3(f_edge.p[1].x - f_edge.p[0].x, f_edge.p[1].y - f_edge.p[0].y, f_edge.p[1].z - f_edge.p[0].z);
					}

					if (external.dotproduct(ab) < 0)
						f_edge.convex = 1;
					else
						f_edge.convex = 1;

					feature_edges.push_back(f_edge);
				}
			}
		}
		else
			printf("CW oriented\n");
	}
}

void find_cutting_path(stl_file *stl, std::vector<feature_edge> &feature_edges)
{
	int i, j, l, m, k, v1, v2, found, start[2], neighbors, neighbors_neighbor, highestacosthetafacet = 0, highestacosthetafacetneigh = 0;
	float acostheta, highestacostheta = 0.0f;
	std::set<int> seen_facets;
	facetstheta *all_facetstheta;
	bool first = true, notclosed;
	feature_edge f_edge, f_edge_previous;
	stl_vertex tmp;

	all_facetstheta = (facetstheta*)malloc(stl->stats.number_of_facets * sizeof(facetstheta));

	for (i = 0; i < stl->stats.number_of_facets; i++) {
		for (j = 0; j < 3; j++) {
			neighbors = stl->neighbors_start[i].neighbor[j];
			for (found = 0, v1 = 0; v1 < 3 && found < 2; v1++) {
				for (v2 = 0; v2 < 3 && found < 2; v2++) {
					if (!memcmp(&stl->facet_start[i].vertex[v1], &stl->facet_start[neighbors].vertex[v2], sizeof(stl_vertex))) {
						f_edge.p[found] = stl->facet_start[i].vertex[v1];
						found++;
						break;
					}

				}
			}
			all_facetstheta[i].start[j] = f_edge.p[0];
			all_facetstheta[i].end[j] = f_edge.p[1];
			acostheta = acos(fabs(mymin(stl->facet_start[i].normal.x*stl->facet_start[neighbors].normal.x +
				stl->facet_start[i].normal.y*stl->facet_start[neighbors].normal.y +
				stl->facet_start[i].normal.z*stl->facet_start[neighbors].normal.z, 1.0f)));
			all_facetstheta[i].neigh[j] = acostheta;
			if (acostheta > highestacostheta) {
				highestacostheta = acostheta;
				highestacosthetafacet = i;
				highestacosthetafacetneigh = j;
			}
		}
	}


	if (highestacostheta < DEG_TO_RAD(10.0)) {
		printf("No cutting path found\n");
		return;
	}


	notclosed = true;
	while (notclosed) {
		// the most critical cosine has been computed, start with lowestcostthetafacet and build a path

		// path may not use these facets again
		seen_facets.insert(highestacosthetafacet);
		seen_facets.insert(stl->neighbors_start[highestacosthetafacet].neighbor[highestacosthetafacetneigh]);

		// find the common edge
		f_edge.facet1_number = highestacosthetafacet;
		f_edge.facet1_neighbor = highestacosthetafacetneigh;
		f_edge.facet2_number = stl->neighbors_start[highestacosthetafacet].neighbor[highestacosthetafacetneigh];
		f_edge.p[0] = all_facetstheta[highestacosthetafacet].start[highestacosthetafacetneigh];
		f_edge.p[1] = all_facetstheta[highestacosthetafacet].end[highestacosthetafacetneigh];
		f_edge.convex = 0;

		if (!first && memcmp(&f_edge.p[0], &f_edge_previous.p[1], sizeof(stl_vertex))) {
			// we want a path..
			tmp = f_edge.p[0];
			f_edge.p[0] = f_edge.p[1];
			f_edge.p[1] = tmp;
		}

		first = false;

		feature_edges.push_back(f_edge);

		f_edge_previous = f_edge;


		// get next following the neighbor with high costheta value
		start[0] = highestacosthetafacet;
		start[1] = stl->neighbors_start[highestacosthetafacet].neighbor[highestacosthetafacetneigh];
		highestacostheta = -1.0f;
		for (l = 0; l < 2; l++) {
			for (j = 0; j < 3; j++) {
				neighbors = stl->neighbors_start[start[l]].neighbor[j];
				if (seen_facets.find(neighbors) != seen_facets.end())continue;
				seen_facets.insert(neighbors);
				for (k = 0; k < 3; k++) {
					neighbors_neighbor = stl->neighbors_start[neighbors].neighbor[k];
					if (seen_facets.find(neighbors_neighbor) != seen_facets.end())continue;
					//seen_facets.insert(neighbors_neighbor);					
					for (m = 0; m < 3; m++) {
						if (memcmp(&f_edge.p[1], &all_facetstheta[neighbors_neighbor].start[m], sizeof(stl_vertex)) &&
							memcmp(&f_edge.p[1], &all_facetstheta[neighbors_neighbor].end[m], sizeof(stl_vertex)))continue;
						if (all_facetstheta[neighbors_neighbor].neigh[m] > highestacostheta) {
							highestacostheta = all_facetstheta[neighbors_neighbor].neigh[m];
							highestacosthetafacet = neighbors_neighbor;
							highestacosthetafacetneigh = m;
							//f_edge_next.facet1_number = next_neigh;
							//f_edge_next.facet2_number = stl->neighbors_start[next_neigh].neighbor[m];
							//f_edge_next.p[0] = f_edge.p[1];
						}
					}
				}
			}
		}

		if (highestacostheta < 0) {// no continuation found, so start a new cutting path
			for (i = 0; i < stl->stats.number_of_facets; i++) {
				if (seen_facets.find(i) != seen_facets.end())continue;
				for (j = 0; j < 3; j++) {
					neighbors = stl->neighbors_start[i].neighbor[j];
					if (all_facetstheta[i].neigh[j] > highestacostheta) {
						highestacostheta = all_facetstheta[i].neigh[j];
						highestacosthetafacet = i;
						highestacosthetafacetneigh = neighbors;
					}
				}
				break;
				//return;
			}
			//notclosed = false;
			// we should push the path and start a new one
			if (i >= stl->stats.number_of_facets || highestacostheta < DEG_TO_RAD(10.0))
				notclosed = false;
		}
	}

	free(all_facetstheta);
}


void find_edge_path(stl_file *stl, std::vector<feature_edge> &feature_edges)
{
	stl_facet facet;
	stl_hash_edge edge;
	stl_hash_edge *link;
	feature_edge f_edge;
	int i, j;

	if (stl->error) return;

	stl->stats.connected_edges = 0;
	stl->stats.connected_facets_1_edge = 0;
	stl->stats.connected_facets_2_edge = 0;
	stl->stats.connected_facets_3_edge = 0;

	my_stl_initialize_facet_check_exact(stl);

	for (i = 0; i < stl->stats.number_of_facets; i++) {
		facet = stl->facet_start[i];
		for (j = 0; j < 3; j++) {
			edge.facet_number = i;
			edge.which_edge = j;
			my_stl_load_edge_exact(stl, &edge, &facet.vertex[j],
				&facet.vertex[(j + 1) % 3]);

			my_insert_hash_edge(stl, edge);
		}
	}

	// we have a list of edges...

	for (i = 0; i < stl->M; i++) {
		link = stl->heads[i];
		while (link != stl->tail) {
			if (link->which_edge > 2) {
				memcpy(&f_edge.p[1], &link->key[0], sizeof(stl_vertex));
				memcpy(&f_edge.p[0], &link->key[3], sizeof(stl_vertex));
			}
			else {
				memcpy(&f_edge.p[0], &link->key[0], sizeof(stl_vertex));
				memcpy(&f_edge.p[1], &link->key[3], sizeof(stl_vertex));
			}
			f_edge.facet1_number = link->facet_number;
			f_edge.facet1_neighbor = link->which_edge;
			f_edge.facet2_number = stl->neighbors_start[link->facet_number].neighbor[link->which_edge % 3];
			f_edge.convex = 1;
			f_edge.acos = acos(mymax(mymin(stl->facet_start[f_edge.facet1_number].normal.x*stl->facet_start[f_edge.facet2_number].normal.x +
				stl->facet_start[f_edge.facet1_number].normal.y*stl->facet_start[f_edge.facet2_number].normal.y +
				stl->facet_start[f_edge.facet1_number].normal.z*stl->facet_start[f_edge.facet2_number].normal.z, 1.0f), -1.0f));
			f_edge.is_ccw = (double)(stl->facet_start[f_edge.facet1_number].vertex[2].x - stl->facet_start[f_edge.facet1_number].vertex[1].x)
				*(double)(stl->facet_start[f_edge.facet1_number].vertex[0].y - stl->facet_start[f_edge.facet1_number].vertex[1].y)
				- (double)(stl->facet_start[f_edge.facet1_number].vertex[2].y - stl->facet_start[f_edge.facet1_number].vertex[1].y)
				*(double)(stl->facet_start[f_edge.facet1_number].vertex[0].x - stl->facet_start[f_edge.facet1_number].vertex[1].x);

			if (f_edge.acos > DEG_TO_RAD(5.0f)) // && acos(f_edge.acos) < DEG_TO_RAD(160))
				feature_edges.push_back(f_edge);
			link = link->next;
		}
	}

	my_stl_free_edges(stl);
}

