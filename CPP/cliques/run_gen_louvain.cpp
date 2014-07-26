#include "stability_gen.h"

#include "louvain_gen.h"
#include "vector_partition.h"
#include "io.h"

#include <lemon/smart_graph.h>
#include <vector>

// Call generic Louvain optimisation from command line
// Inputs: graph, null vectors (both as txt files), resolution parameter (time) and random seed (both optional)
// Example: ./run_gen_louvain.sh graph_file null_model_file time random_seed
//
// TODO create run_gen_louvain_ntimes for multiple runs at once
// create specialised version when no null model is given / can be extracted from matrix!?
int main(int argc, char *argv []) {

    // initialise random seed from input
    if (argc < 5) {
        srand(time(0));
    } else {
        srand(atoi(argv[4]));
    }

    // initialise markov time from input
    double current_markov_time = 1;
    if (argc > 3 ){
        current_markov_time = atof(argv[3]);
    }

    // initialise graph and edge map
    lemon::SmartGraph input_graph;
    lemon::SmartGraph::EdgeMap<double> input_graph_weights(input_graph);

    
    clq::read_edgelist_weighted_graph(argv[1], input_graph, input_graph_weights);
    
    int num_nodes = lemon::countNodes(input_graph);
    std::vector<std::vector<double>> null_model = clq::read_null_model(argv[2],num_nodes);
    //clq::print_2d_vector(null_model);

    
    // define type for vector partition initialise start partition, qualitiy functions etc.
    typedef clq::VectorPartition partition;

    partition start_partition(lemon::countNodes(input_graph));
    start_partition.initialise_as_singletons();

    clq::find_linearised_generic_stability quality(current_markov_time);
    clq::linearised_generic_stability_gain quality_gain(current_markov_time);
    
    std::vector<partition> optimal_partitions;

    // call Louvain
    clq::output("Start Louvain");
    double stability = clq::find_optimal_partition_louvain_gen<partition>(
                           input_graph, input_graph_weights, null_model, quality, quality_gain,
                           start_partition, optimal_partitions, 1e-12);
    
    // store output partitions in file
    clq::partitions_to_file("optimal_partitions.dat", optimal_partitions);

    // display final result
    partition best_partition = optimal_partitions.back();
    clq::print_partition_list(best_partition);
    clq::output("Stability: ", stability);

    return 0;
}
