#pragma once

#include <vector>
#include <lemon/concepts/graph.h>

#include "graphhelpers.h"
#include "linearised_internals_generic.h"
#include "generators.h"

// TODO: Separate out louvain into smaller functions


namespace clq {

/**
 @brief  Louvain method - greedy algorithm to find community structure of a network.

 The second phase of the algorithm consists in building a new network
 whose nodes are now the communities found during the first phase.
 To do so, the weights of the links between the new nodes are given by
 the sum of the weight of the links between nodes in the corresponding
 two communities.
 Links between nodes of the same community lead to self-loops for this
 community in the new network. Once this second phase is completed,
 it is then possible to reapply the first phase of the algorithm
 to the resulting weighted network and to iterate.

 @param[in]  graph     graph to partition
 @param[in]  weights   EdgeMap of edge weights
 @param[in]  compute_quality     partition quality functor
 @param[out] optimal_partitions     optimal partitions found, last in vector is overall best
 */
template<typename P, typename T, typename W, typename QF, typename QFDIFF>
double find_optimal_partition_louvain_gen(T &graph, W &weights,
        std::vector<std::vector<double>> null_model_vec, QF compute_quality,
        QFDIFF compute_quality_diff, P initial_partition,
        std::vector<P> &optimal_partitions, double minimum_improve) {

    typedef typename T::Node Node;
    typedef typename T::Edge Edge;
    typedef typename T::NodeIt NodeIt;
    typedef typename T::EdgeIt EdgeIt;
    typedef typename T::IncEdgeIt IncEdgeIt;

    P partition(lemon::countNodes(graph));
    partition = initial_partition;
    P partition_init = initial_partition;
    if (!optimal_partitions.empty()) {
        partition_init = optimal_partitions.back();
    }
    auto internals = clq::gen(compute_quality, graph, weights, partition, partition_init, null_model_vec);

    //double minimum_improve = 0.000000001; //1e-9
    double current_quality = compute_quality(internals);
    //clq::output("current_quality", current_quality);
    bool one_level_end = false;
    double old_quality = current_quality;
    bool did_nodes_move = false;

    // Randomise the looping over nodes
    // initialise random number generator outside louvain when calling externally multiple times!
    // srand(std::time(0));
    // create vector to shuffle
    std::vector<Node> nodes_ordered_randomly;
    for (NodeIt temp_node(graph); temp_node != lemon::INVALID; ++temp_node) {
        nodes_ordered_randomly.push_back(temp_node);
    }
    // Randomly shuffle nodes
    std::random_shuffle(nodes_ordered_randomly.begin(),
                        nodes_ordered_randomly.end());

    do {
        did_nodes_move = false;
        old_quality = current_quality;

        // loop over all nodes in random order
        typename std::vector<Node>::iterator n1_it;
        for (n1_it = nodes_ordered_randomly.begin(); n1_it
                != nodes_ordered_randomly.end(); ++n1_it) {

            // get node id and comm id
            Node n1 = *n1_it;
            unsigned int node_id = graph.id(n1);
            unsigned int comm_id = partition.find_set(node_id);
            isolate_and_update_internals(graph, weights, n1, internals,
                                         partition);
            //default option for re-inclusion of node
            unsigned int best_comm = comm_id;
            double best_gain = 0;

            unsigned int num_neighbour_communities =
                internals.neighbouring_communities_list.size();

            // loop over neighbouring communities
            for (unsigned int k = 0; k < num_neighbour_communities; ++k) {

                unsigned int comm_id_neighbour =
                    internals.neighbouring_communities_list[k];
                double gain = compute_quality_diff(internals,
                                                   comm_id_neighbour, node_id);
                if (gain > best_gain) {
                    best_comm = comm_id_neighbour;
                    best_gain = gain;
                    // avoid not necessary movements, place node in old community if possible
                } else if (gain == best_gain && comm_id == comm_id_neighbour) {
                    best_comm = comm_id;
                }

            }

            insert_and_update_internals(graph, weights, n1, internals,
                                        partition, best_comm);


            if (best_comm != comm_id) {
                did_nodes_move = true;
            }
        }
        if (did_nodes_move == true) {
            current_quality = compute_quality(internals);
            one_level_end = true;
        }

        // If there was a movement of the nodes AND quality increases make another run
    } while ((current_quality - old_quality) > minimum_improve);

    //clq::output("current_quality", current_quality);
    // Start Second phase - create reduced graph with self loops
    std::map<int, int> new_comm_id_to_old_comm_id = partition.normalise_ids();

    int hierarchy = optimal_partitions.size();
    if (one_level_end == true) {
        // Compile P into original partition size
        if (hierarchy == 0) {
            optimal_partitions.push_back(partition);
        } else {
            // get size of partition one level below, i.e. number of nodes in original graph
            unsigned int original_number_of_nodes =
                optimal_partitions[hierarchy - 1].element_count();
            // create new empty partition of this size
            P partition_original_nodes(original_number_of_nodes);

            // loop over nodes one level below
            int old_comm, new_comm;
            for (unsigned int id = 0; id < original_number_of_nodes; id++) {
                // get the communities for each node one level below
                old_comm = optimal_partitions[hierarchy - 1].find_set(id);
                // use this as node_id in the current partition as old community id
                // is equivalent to new node id and read out new community id
                new_comm = partition.find_set(old_comm);
                // include pair (node, new community) id in the newly created partition
                partition_original_nodes.add_node_to_set(id, new_comm);
            }
            optimal_partitions.push_back(partition_original_nodes);
        }

        // Create graph from partition
        T reduced_graph;
        W reduced_weight_map(reduced_graph);
        create_reduced_graph_from_partition(reduced_graph, reduced_weight_map,
                                            graph, weights, partition, new_comm_id_to_old_comm_id,
                                            internals);

        // get number of nodes in new reduced graph and initialise partition + new null model vectors
        unsigned int num_nodes_reduced_graph = lemon::countNodes(reduced_graph);
        P singleton_partition(num_nodes_reduced_graph);
        singleton_partition.initialise_as_singletons();

        // initialise reduced null model_vec
        std::vector<std::vector<double>> reduced_null_model_vec( null_model_vec.size(), std::vector<double>(num_nodes_reduced_graph,0) );

        for (unsigned int k = 0; k<num_nodes_reduced_graph; ++k)
        {
            for (unsigned int j =0; j< null_model_vec.size(); ++j) 
            {   
                // reduced model null vector correspond to summed loss terms per communitiy in the old graph..
                reduced_null_model_vec[j][k] = internals.comm_loss_vectors[j][new_comm_id_to_old_comm_id[k]];
            }
        }
        //clq::output("reduced null vectors");
        //clq::print_2d_vector(reduced_null_model_vec);

        return clq::find_optimal_partition_louvain_gen<P>(reduced_graph,
                reduced_weight_map, reduced_null_model_vec, compute_quality,
                compute_quality_diff, singleton_partition, optimal_partitions,
                minimum_improve);
    }
    if (hierarchy == 0) {
        optimal_partitions.push_back(partition);
    }
    return current_quality;
}

}
