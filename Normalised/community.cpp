#include "community.h"
extern "C" {
#include <math.h>
#include "mex.h"
}

using namespace std;

Community::~Community() {
	free(g.links);
	free(g.weights);
	free(g.degrees);
}

Community::Community(char * filename, char * filename_w, int type, int nbp,
		double minm) {
	g = Graph(filename, filename_w, type);
	size = g.nb_nodes;

	neigh_weight.resize(size, -1);
	neigh_pos.resize(size);
	neigh_last = 0;

	n2c.resize(size);
	in.resize(size);
	tot.resize(size);

	for (int i = 0; i < size; i++) {
		n2c[i] = i;
		in[i] = g.nb_selfloops(i);
		tot[i] = g.weighted_degree(i);
	}

	nb_pass = nbp;
	min_modularity = minm;

}

Community::Community(double * data, int length_data, int nbp, double minm,
		double timet, int type) {

	g = Graph(data, length_data, type);
	size = g.nb_nodes;

	neigh_weight.resize(size, -1);
	neigh_pos.resize(size);
	neigh_last = 0;

	n2c.resize(size);
	in.resize(size);
	tot.resize(size);

	for (int i = 0; i < size; i++) {
		n2c[i] = i;
		in[i] = g.nb_selfloops(i);
		tot[i] = g.weighted_degree(i);
	}

	nb_pass = nbp;
	min_modularity = minm;

	time = timet;

}

Community::Community(Graph gc, int nbp, double minm, double timet) {
	g = gc;
	size = g.nb_nodes;

	//cout << "SIZE = " << size << endl;
	neigh_weight.resize(size, -1);
	neigh_pos.resize(size);
	neigh_last = 0;

	n2c.resize(size);
	in.resize(size);
	tot.resize(size);

	for (int i = 0; i < size; i++) {
		n2c[i] = i;
		in[i] = g.nb_selfloops(i);
		tot[i] = g.weighted_degree(i);
	}

	nb_pass = nbp;
	min_modularity = minm;

	time = timet;
}

void Community::display() {
	cerr << endl << "<";
	for (int i = 0; i < size; i++)
		cerr << " " << i << "/" << n2c[i] << "/" << in[i] << "/" << tot[i];
	cerr << ">" << endl;
}

double Community::modularity() {
	double q = 1.0 - time;
	double m2 = (double) g.total_weight;

	for (int i = 0; i < size; i++) {
		if (tot[i] > 0)
			q += time * (double) in[i] / m2 - ((double) tot[i] / m2)
					* ((double) tot[i] / m2);
	}

	return q;
}

void Community::neigh_comm(unsigned int node) {
	for (unsigned int i = 0; i < neigh_last; i++)
		neigh_weight[neigh_pos[i]] = -1;
	neigh_last = 0;

	pair<unsigned int *, double *> p = g.neighbors(node);

	unsigned int deg = g.nb_neighbors(node);

	neigh_pos[0] = n2c[node];
	neigh_weight[neigh_pos[0]] = 0;
	neigh_last = 1;

	for (unsigned int i = 0; i < deg; i++) {
		unsigned int neigh = *(p.first + i);
		unsigned int neigh_comm = n2c[neigh];
		double neigh_w = (g.weights == NULL) ? 1. : *(p.second + i);

		if (neigh != node) {
			if (neigh_weight[neigh_comm] == -1) {
				neigh_weight[neigh_comm] = 0.;
				neigh_pos[neigh_last++] = neigh_comm;
			}
			neigh_weight[neigh_comm] += neigh_w;
		}
	}
}

void Community::partition2graph() {
	vector<int> renumber(size, -1);
	for (int node = 0; node < size; node++) {
		renumber[n2c[node]]++;
	}

	int final = 0;
	for (int i = 0; i < size; i++)
		if (renumber[i] != -1)
			renumber[i] = final++;

	for (int i = 0; i < size; i++) {
		pair<unsigned int *, double *> p = g.neighbors(i);

		int deg = g.nb_neighbors(i);
		for (int j = 0; j < deg; j++) {
			int neigh = *(p.first + j);
			cout << renumber[n2c[i]] << " " << renumber[n2c[neigh]] << endl;
		}
	}
}

void Community::display_partition(char* outfilename) {
	printf("displayin");
	ofstream foutput(outfilename, ios::out | ios::app);

	vector<int> renumber(size, -1);
	for (int node = 0; node < size; node++) {
		renumber[n2c[node]]++;
	}

	int final = 0;
	for (int i = 0; i < size; i++)
		if (renumber[i] != -1)
			renumber[i] = final++;

	for (int i = 0; i < size; i++) {
		foutput << i << " " << renumber[n2c[i]] << endl;
		printf("%d", i);
		printf(" ");
		printf("%d", renumber[n2c[i]]);
		printf("\n");
	}
	printf("displayout");
}

vector<vector<int> > Community::display_partition2(vector<vector<int> > output) {

	// renumber nodes
	// first mark communities
	vector<int> renumber(size, -1);
	for (int node = 0; node < size; node++) {
		renumber[n2c[node]]++;
	}
	// then renumber communities
	int final = 0;
	for (int i = 0; i < size; i++)
		if (renumber[i] != -1)
			renumber[i] = final++;

	vector<int> temp(size, 0);
	// assign new labels
	for (int i = 0; i < size; i++) {
		temp[i] = renumber[n2c[i]];
	}
	// push to vector
	output.push_back(temp);

	return output;
}

// This function is not so nice
// malloc is dirty
Graph Community::partition2graph_binary() {

	// Renumber communities
	vector<int> renumber(size, -1);
	for (int node = 0; node < size; node++) {
		renumber[n2c[node]]++;
	}

	int final = 0;
	for (int i = 0; i < size; i++)
		if (renumber[i] != -1)
			renumber[i] = final++;

	// Compute communities
	vector<vector<int> > comm_nodes(final);
	for (int node = 0; node < size; node++) {
		comm_nodes[renumber[n2c[node]]].push_back(node);
	}

	// Compute weighted graph
	Graph g2;
	g2.nb_nodes = comm_nodes.size();
	g2.degrees = (unsigned long *) malloc(comm_nodes.size() * sizeof(long));
	g2.links = (unsigned int *) malloc((long) g.nb_links * sizeof(int));
	g2.weights = (double *) malloc((long) g.nb_links * sizeof(double));

	long where = 0;
	int comm_deg = comm_nodes.size();
	for (int comm = 0; comm < comm_deg; comm++) {
		map<int, double> m;
		map<int, double>::iterator it;

		int comm_size = comm_nodes[comm].size();
		for (int node = 0; node < comm_size; node++) {
			pair<unsigned int *, double *> p = g.neighbors(
					comm_nodes[comm][node]);
			int deg = g.nb_neighbors(comm_nodes[comm][node]);
			for (int i = 0; i < deg; i++) {
				int neigh = *(p.first + i);
				int neigh_comm = renumber[n2c[neigh]];
				double neigh_weight = (g.weights == NULL) ? 1.
						: *(p.second + i);

				it = m.find(neigh_comm);
				if (it == m.end())
					m.insert(make_pair(neigh_comm, neigh_weight));
				else
					it->second += neigh_weight;
			}
		}

		g2.degrees[comm] = (comm == 0) ? m.size() : g2.degrees[comm - 1]
				+ m.size();
		g2.nb_links += m.size();

		for (it = m.begin(); it != m.end(); it++) {
			g2.total_weight += it->second;
			g2.links[where] = it->first;
			g2.weights[where] = it->second;
			where++;
		}
		//    cout << comm << " " << g2.weighted_degrees[comm] << endl;
	}

	g2.links = (unsigned int*) realloc(g2.links, (long) g2.nb_links
			* sizeof(int));
	g2.weights = (double*) realloc(g2.weights, (long) g2.nb_links
			* sizeof(double));

	return g2;
}


bool Community::one_level() {
	bool improvement = false;
	int nb_moves;
	int nb_pass_done = 0;
	double new_mod = modularity();
	double cur_mod = new_mod;

	// repeat while
	//   there is an improvement of modularity
	//   or there is an improvement of modularity greater than a given epsilon
	//   or a predefined number of pass have been done


	vector<int> random_order(size);
	for (int i = 0; i < size; i++)
		random_order[i] = i;
	for (int i = 0; i < size - 1; i++) {
		int rand_pos = rand() % (size - i) + i;
		int tmp = random_order[i];
		random_order[i] = random_order[rand_pos];
		random_order[rand_pos] = tmp;
	}

	do {
		cur_mod = new_mod;
		nb_moves = 0;
		nb_pass_done++;

		// for each node: remove the node from its community and insert it in the best community
		for (int node_tmp = 0; node_tmp < size; node_tmp++) {
			int node = random_order[node_tmp];
			int node_comm = n2c[node];
			double w_degree = g.weighted_degree(node);

			// computation of all neighboring communities of current node
			neigh_comm(node);
			// remove node from its current community
			remove(node, node_comm, neigh_weight[node_comm]);

			// compute the nearest community for node
			// default choice for future insertion is the former community
			int best_comm = node_comm;
			double best_nblinks = 0.;
			double best_increase = 0.;
			for (unsigned int i = 0; i < neigh_last; i++) {
				double increase = modularity_gain(node, neigh_pos[i],
						neigh_weight[neigh_pos[i]], w_degree);
				if (increase > best_increase) {
					best_comm = neigh_pos[i];
					best_nblinks = neigh_weight[neigh_pos[i]];
					best_increase = increase;
				}
			}

			// insert node in the nearest community
			insert(node, best_comm, best_nblinks);

			if (best_comm != node_comm)
				nb_moves++;

		}

		new_mod = modularity();

		if (nb_moves > 0)
			improvement = true;

		if (new_mod > 1) {
			printf(
					"\tThe louvain algorithm appears to be stuck in a loop. \n\tPlease first make sure that the adjacency matrix is perfectly symmetric \n\tand, if it is, try to increase the precision.\n");
			for (int i = 0; i < size; i++)
				n2c[i] = 0;
			nb_pass = -10;
			return false;
		}

	} while (nb_moves > 0 && new_mod - cur_mod > min_modularity);

	return improvement;

}

