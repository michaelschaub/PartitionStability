

#if defined(macintosh) || defined(MACOS) || defined(_MACOS)  || defined(__APPLE__)  || defined(_MAC)   || defined(MAC)  || defined(mac)  || defined(MACINTOSH)
#define __MAC__
#endif

#if defined(_WIN32)  || defined(WIN32)  || defined(__WIN32__) || defined(_WIN64)  || defined(WIN64) || defined(__WIN64__) || defined(_WINDOWS) || defined(__WINDOWS__)
#define __WIN__
#endif

#if defined(linux)  || defined(__linux__)  || defined(_linux) || defined(LINUX)  || defined(_LINUX) || defined(_UNIX) || defined(__UNIX__) || defined(__gnu_linux__) || defined(__unix__) || defined(UNIX) || defined(unix) || defined(sparc)
#define __lin__
#endif

#ifndef GRAPH_H
#define GRAPH_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>

#ifdef __lin__
#include <malloc.h>
#endif

#ifdef __WIN__
#include <malloc.h>
#endif

#ifdef __MAC__
#include <malloc/malloc.h>
#endif

#define WEIGHTED   0
#define UNWEIGHTED 1

using namespace std;

class Graph {
public:
	unsigned int nb_nodes;
	unsigned int nb_links;
	double total_weight;

	unsigned long *degrees;
	unsigned int *links;
	double *weights;

	Graph();

	// binary file format is
	// 4 bytes for the number of nodes in the graph
	// 4*(nb_nodes) bytes for the cumulative degree for each node:
	//    deg(0)=degrees[0]
	//    deg(k)=degrees[k]-degrees[k-1]
	// 4*(sum_degrees) bytes for the links
	// IF WEIGHTED 4*(sum_degrees) bytes for the weights
	Graph(char *filename, char *filename_w, int type);

	Graph(double * data, int length_data, int type);

	Graph(int nb_nodes, int nb_links, double total_weight, int *degrees,
			int *links, double *weights);

	void display(void);
	void display_reverse(void);
	void display_binary(char *outfile);
	bool check_symmetry();

	// return the number of neighbors (degree) of the node
	inline unsigned int nb_neighbors(unsigned int node);

	// return the number of self loops of the node
	inline double nb_selfloops(unsigned int node);

	// return the weighted degree of the node
	inline double weighted_degree(unsigned int node);

	// return pointers to the first neighbor and first weight of the node
	inline pair<unsigned int *, double *> neighbors(unsigned int node);
};

inline unsigned int Graph::nb_neighbors(unsigned int node) {
	assert(node>=0 && node<nb_nodes);

	if (node == 0)
		return degrees[0];
	else
		return degrees[node] - degrees[node - 1];
}

inline double Graph::nb_selfloops(unsigned int node) {
	assert(node>=0 && node<nb_nodes);

	pair<unsigned int *, double *> p = neighbors(node);
	for (unsigned int i = 0; i < nb_neighbors(node); i++) {
		if (*(p.first + i) == node) {
			if (weights != NULL)
				return (double) *(p.second + i);
			else
				return 1.;
		}
	}
	return 0;
}

inline double Graph::weighted_degree(unsigned int node) {
	assert(node>=0 && node<nb_nodes);

	if (weights == NULL)
		return (double) nb_neighbors(node);
	else {
		pair<unsigned int *, double *> p = neighbors(node);
		double res = 0;
		for (unsigned int i = 0; i < nb_neighbors(node); i++) {
			res += (double) *(p.second + i);
		}
		return res;
	}
}

inline pair<unsigned int *, double *> Graph::neighbors(unsigned int node) {
	assert(node>=0 && node<nb_nodes);

	if (node == 0)
		return make_pair(links, weights);
	else if (weights != NULL)
		return make_pair(links + (long) degrees[node - 1], weights
				+ (long) degrees[node - 1]);
	else
		return make_pair(links + (long) degrees[node - 1], weights);
}

#endif // GRAPH_H
