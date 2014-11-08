#include <fstream>
#include "graph_binary.h"
#include "math.h"
extern "C" {
#include <math.h>
#include "mex.h"
#include "matrix.h"
}

#if defined(macintosh) || defined(MACOS) || defined(_MACOS)  || defined(__APPLE__)  || defined(_MAC)   || defined(MAC)  || defined(mac)  || defined(MACINTOSH)
#define __MAC__
#endif

#if defined(_WIN32)  || defined(WIN32)  || defined(__WIN32__) || defined(_WIN64)  || defined(WIN64) || defined(__WIN64__) || defined(_WINDOWS) || defined(__WINDOWS__)
#define __WIN__
#endif

#if defined(linux)  || defined(__linux__)  || defined(_linux) || defined(LINUX)  || defined(_LINUX) || defined(_UNIX) || defined(__UNIX__) || defined(__gnu_linux__) || defined(__unix__) || defined(UNIX) || defined(unix) || defined(sparc)
#define __lin__
#endif

#ifdef __WIN__
#endif

#ifdef __lin__
#include <sys/mman.h>
#endif

#ifdef __MAC__
#include <sys/mman.h>
#endif

Graph::Graph() {
	nb_nodes = 0;
	nb_links = 0;
	total_weight = 0;
}

Graph::Graph(char *filename, int type) {
	ifstream finput;
	finput.open(filename, fstream::in | fstream::binary);

	// read number of nodes on 4 bytes
	finput.read((char *) &nb_nodes, 4);

	degrees = (unsigned long *) malloc((long) nb_nodes * sizeof(long));
	finput.read((char *) degrees, (long) nb_nodes * sizeof(int));

	// read links: 4 bytes for each link (each link is counted twice)
	nb_links = degrees[nb_nodes - 1] / 2;
	links = (unsigned int *) malloc((long) nb_links * 2 * sizeof(int));
	finput.read((char *) links, (long) nb_links * 2 * sizeof(int));

	// IF WEIGHTED : read weights: 4 bytes for each link (each link is counted twice)
	if (type == WEIGHTED) {
		weights = (double *) malloc((long) nb_links * 8);
		finput.read((char *) weights, (long) nb_links * 8);
		total_weight = 0;
		for (unsigned int i = 0; i < nb_links * 2; i++) {
			total_weight += weights[i];
		}
	} else {
		weights = NULL;
		total_weight = 2 * nb_links;
	}

	// New attribute being the backup of the number of nodes per community
	nb_nodes_per_comm = (int *) malloc((long) nb_nodes * 4);
	for (unsigned int i = 0; i < nb_nodes; i++)
		nb_nodes_per_comm[i] = 1;

}

Graph::Graph(double * data, int length_data, int type) {

	nb_nodes = int(data[length_data - 1]) + 1;

	degrees = (unsigned long *) malloc((long) nb_nodes * sizeof(unsigned long));

	int tmp_node = -1;
	int j = 0;
	int tot = 0;
	for (int i = 0; i < length_data; i++) {
		tot += 1;
		if (int(data[i]) == tmp_node)
			degrees[tmp_node] += 1;
		else {
			tmp_node++;
			degrees[tmp_node] = tot;
		}
	}

	nb_links = degrees[nb_nodes - 1] / 2;

	links = (unsigned int *) malloc((long) nb_links * 2 * sizeof(unsigned int));
	for (int i = 0; i < length_data; i++) {
		links[i] = int(data[length_data + i]);
	}

	// IF WEIGHTED : read weights: 4 bytes for each link (each link is counted twice)
	weights = NULL;
	total_weight = length_data;
	if (type == WEIGHTED) {
		weights = (double *) malloc((long) nb_links * 2 * sizeof(double));//new double[nb_links];
		total_weight = 0;
		for (unsigned int i = 0; i < nb_links * 2; i++) {
			weights[i] = (double) data[2 * length_data + i];
			total_weight = total_weight + weights[i];
		}
	}

	// New attribute being the backup of the number of nodes per community
	nb_nodes_per_comm = (int *) malloc((long) nb_nodes * sizeof(int));
	for (unsigned int i = 0; i < nb_nodes; i++)
		nb_nodes_per_comm[i] = 1;

}

void Graph::display() {
	for (unsigned int node = 0; node < nb_nodes; node++) {
		pair<unsigned int *, double *> p = neighbors(node);
		for (unsigned int i = 0; i < nb_neighbors(node); i++) {
			if (node <= *(p.first + i)) {
				if (weights != NULL)
					cout << node << " " << *(p.first + i) << " " << *(p.second
							+ i) << endl;
				else
					cout << node << " " << *(p.first + i) << endl;
			}
		}
	}
}

void Graph::display_reverse() {
	for (unsigned int node = 0; node < nb_nodes; node++) {
		pair<unsigned int *, double *> p = neighbors(node);
		for (unsigned int i = 0; i < nb_neighbors(node); i++) {
			if (node > *(p.first + i)) {
				if (weights != NULL)
					cout << *(p.first + i) << " " << node << " " << *(p.second
							+ i) << endl;
				else
					cout << *(p.first + i) << " " << node << endl;
			}
		}
	}
}

bool Graph::check_symmetry() {
	int error = 0;
	for (unsigned int node = 0; node < nb_nodes; node++) {
		pair<unsigned int *, double *> p = neighbors(node);
		for (unsigned int i = 0; i < nb_neighbors(node); i++) {
			unsigned int neigh = *(p.first + i);
			double weight = *(p.second + i);

			pair<unsigned int *, double *> p_neigh = neighbors(neigh);
			for (unsigned int j = 0; j < nb_neighbors(neigh); j++) {
				unsigned int neigh_neigh = *(p_neigh.first + j);
				double neigh_weight = *(p_neigh.second + j);

				if (node == neigh_neigh && weight != neigh_weight) {
					cout << node << " " << neigh << " " << weight << " "
							<< neigh_weight << endl;
					if (error++ == 10)
						exit(0);
				}
			}
		}
	}
	return (error == 0);
}

void Graph::display_binary(char *outfile) {
	ofstream foutput;
	foutput.open(outfile, fstream::out | fstream::binary);

	foutput.write((char *) (&nb_nodes), 4);
	foutput.write((char *) (degrees), 4 * nb_nodes);
	foutput.write((char *) (links), 8 * nb_links);
}
