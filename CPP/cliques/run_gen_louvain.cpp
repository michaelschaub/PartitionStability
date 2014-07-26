#include "stability_gen.h"

#include "louvain_gen.h"
#include "vector_partition.h"
#include "io.h"

#include <lemon/smart_graph.h>
#include <vector>

// Call generic Louvain optimisation from command line
// Inputs: graph, null vectors, resolution parameter
int main(int argc, char *argv []) {
    // define type for vector partition
    typedef clq::VectorPartition partition;
    lemon::SmartGraph input_graph;
    lemon::SmartGraph::EdgeMap<double> input_graph_weights(input_graph);

    double stability = 0;
    double current_markov_time = 1;

    clq::read_edgelist_weighted_graph(argv[1], input_graph, input_graph_weights);
    int num_nodes = lemon::countNodes(input_graph);
    
        
    std::vector<std::vector<double>> null_model =clq::read_null_model(argv[2],num_nodes);
    //null_model[0] = {2.0, 2.0, 3.0, 3.0 ,2.0, 2.0};
    //null_model[1] = {2.0/14, 2.0/14, 3.0/14, 3.0/14 ,2.0/14, 2.0/14};
    clq::print_2d_vector(null_model);

    partition start_partition(lemon::countNodes(input_graph));
    start_partition.initialise_as_singletons();

    clq::find_linearised_generic_stability quality(current_markov_time);
    clq::linearised_generic_stability_gain quality_gain(current_markov_time);

    std::vector<partition> optimal_partitions;

    clq::output("Start Louvain");

    stability = clq::find_optimal_partition_louvain_gen<partition>(
                    input_graph, input_graph_weights, null_model, quality, quality_gain,
                    start_partition, optimal_partitions, 1e-12);
    clq::partitions_to_file("optimal_partitions.dat", optimal_partitions);


    partition best_partition = optimal_partitions.back();
    clq::print_partition_list(best_partition);
    clq::output("Stability: ", stability);

    return 0;
}
