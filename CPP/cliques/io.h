#pragma once
#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <set>
#include <lemon/core.h>
#include <vector>
#include <map>

namespace clq {


//============================================================================================
// PRINT_MATRIX
// Template to print out linear (array/vector) in "matrix" format assuming column first 
// ordering, i.e., A(i,j) corresponds to matrix[i+j*N], where N is the dimension of the matrix
// 
// INPUTS:  matrix -- linear container (array/vector)
//          lda -- dimension of the matrix
//============================================================================================
template <typename T>
void print_matrix(T matrix, int lda) {
    for (int i=0; i<lda; ++i) {
        for (int j=0; j<lda; ++j) {
            std::cout << matrix[j+i*lda] << "\t";
        }
        std::cout << std::endl;
    } 
    std::cout << std::endl;
}



//TODO WRITE DESCRIPTION
std::vector<double> read_edgelist_weighted(std::string filename) {
    // initialise input stream and strings for readout
    std::ifstream my_file(filename.c_str());
    std::string line;
    std::string mystring;

    // check if file is open
    if (!my_file.is_open()) {
        std::cout << "couldn't open file:" << filename << std::endl;
        exit(1);
    }

    // keep track of size of graph and create internal adjacency list
    // NB: Node numbering starts with 0!
    int max_node_id_seen = -1;    
    std::vector<int> from, to;
    std::vector<double> weight;

    //readout contents from my_file into string, line by line
    while (std::getline(my_file, line)) {

        std::stringstream lineStream(line);
        //readout node id and weights
        std::getline(lineStream, mystring, ' ');
        int node1_id = atoi(mystring.c_str());
        from.push_back(node1_id);
        
        std::getline(lineStream, mystring, ' ');
        int node2_id = atoi(mystring.c_str());
        to.push_back(node2_id);

        std::getline(lineStream, mystring, ' '); 
        double edge_weight = atof(mystring.c_str());
        weight.push_back(edge_weight);

        if (node1_id > max_node_id_seen) {
            max_node_id_seen = node1_id;
        }

        if (node2_id > max_node_id_seen) {
            max_node_id_seen = node2_id;
        }
    }
    // don't forget to close file after readout...
    my_file.close();
    
    // now write adjecency matrix in vector form (row first ordering)
    std::vector<double> Adj((max_node_id_seen+1)*(max_node_id_seen+1),0);
    for (unsigned int i =0; i<to.size(); i++){
        int index = from[i]*(max_node_id_seen+1)+to[i];
        Adj[index] = weight[i];
    }
    
    return Adj;
}

//TODO WRITE DESCRIPTION
template<typename G, typename E>
bool read_edgelist_weighted_graph(std::string filename, G &graph, E &weights) {
    // initialise input stream and strings for readout
    std::ifstream my_file(filename.c_str());
    std::string line;
    std::string mystring;

    // check if file is open
    if (!my_file.is_open()) {
        std::cout << "couldn't open file:" << filename << std::endl;
        exit(1);
    }

    // define Node class for convenience
    typedef typename G::Node Node;
    
    // keep track of graph size
    // NB: numbering starts with 0!
    int max_node_id_seen = -1;

    //readout contents from my_file into string, line by line
    while (std::getline(my_file, line)) {

        std::stringstream lineStream(line);
        //readout node id and weights
        std::getline(lineStream, mystring, ' ');
        int node1_id = atoi(mystring.c_str());
        std::getline(lineStream, mystring, ' ');
        int node2_id = atoi(mystring.c_str());
        std::getline(lineStream, mystring, ' ');
        double weight = atof(mystring.c_str());

        if (node1_id > max_node_id_seen) {
            int difference = node1_id - max_node_id_seen;
            for (int i=0;i<difference;++i) {
                graph.addNode();
            }
            max_node_id_seen = node1_id;
        }

        if (node2_id > max_node_id_seen) {
            int difference = node2_id - max_node_id_seen;
            for (int i=0;i<difference;++i) {
                graph.addNode();
            }
            max_node_id_seen = node2_id;
        }

        Node node1 = graph.nodeFromId(node1_id);
        Node node2 = graph.nodeFromId(node2_id);

        typename G::Edge edge = graph.addEdge(node1, node2);
        weights.set(edge, weight);
    }

    my_file.close();
    return true;
}


//TODO WRITE DESCRIPTION
template<typename G, typename E>
void write_edgelist_weighted(std::string filename, G &graph, E &weights) {
    // initialise input stream and strings for readout
    std::ofstream my_file(filename.c_str()); 
    std::string mystring;

    // check if file is open
    if (!my_file.is_open()) {
        std::cout << "couldn't open file" << std::endl;
        exit(1);
    }

    for(typename G::EdgeIt e(graph); e!=lemon::INVALID; ++e){
        int node1 = graph.id(graph.u(e));
        int node2 = graph.id(graph.v(e));
        double weigth = weights[e];

        my_file << node1 << " " << node2 << " " << weigth << "\n";
    }

    my_file.close();
}

}
