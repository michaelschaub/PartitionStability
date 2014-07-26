#pragma once

#include <vector>
#include <lemon/smart_graph.h>
#include <time.h>
#include "graphhelpers.h"
#include "linearised_internals_generic.h"

namespace clq
{
/**
 @brief  Functor for finding generic stability of weighted graph
 */
struct find_linearised_generic_stability
{
    // Markov time state variable not used so far --- maybe useful in the future, now set to 1 by default.
    double markov_time;

    // constructor
    find_linearised_generic_stability(double markov_time = 1.0) :
        markov_time(markov_time)
    {
    }

    template<typename G, typename P, typename W>
    double operator ()(G &graph, P &partition, W &weights,P &partition_init, std::vector< std::vector<double> > null_model_vectors)
    {
        clq::LinearisedInternalsGeneric internals(graph, weights, partition, partition_init, null_model_vectors);
        return (*this)(internals);
    }

    template<typename I>
    double operator ()(I &internals)
    {
        // first part equals 1-t
        double q = 1-markov_time;
        double q2 = 0;

        // loop over all communities and sum up contributions (gain - loss)
        unsigned int num_null_model_vec = internals.num_null_model_vectors;
        unsigned int num_nodes = internals.num_nodes;
        // loop over all possible community indices
        for (unsigned int i = 0; i < num_nodes; i++)
        {   // for each possible index check if the community is non-empty, i.e. if there is a loss term 
            bool check = 0;
            //  for this to work null model vectors should be passed as pairs!
            for (unsigned int j = 0; j < num_null_model_vec; j=j+2)
            {
                // if there is a contribution sum up negative terms first
                if(internals.comm_loss_vectors[j][i] !=0)
                {
                    check =1;
                    //clq::output("check", check, "\n");
                    //clq::print_2d_vector(internals.comm_loss_vectors);
                    q2 -= double(internals.comm_loss_vectors[j][i])*double(internals.comm_loss_vectors[j+1][i]);
                }
            }

            // if there was a contribution from this community get the gain term as well 
            if (check)
            {
                q2 += markov_time * double(internals.comm_w_in[i]);
            }
        }
 
        return q + q2;
    }

};

/**
 @brief  Functor for finding stability gain (normalised Laplacian) with for weighted graph
 */
struct linearised_generic_stability_gain
{
    //state variable unused so far -- for future use
    double markov_time;

    //constructor
    linearised_generic_stability_gain(double mtime = 1.0) :
        markov_time(mtime)
    {
    }

    template<typename I>
    double operator ()(I &internals, int comm_id_neighbour, int node_id)
    {
        // compute additonal loss factor...
        double comm_loss = 0;
        //clq::output(comm_loss,internals.num_null_model_vectors,"loss\n");
        for (unsigned int i = 0; i < internals.num_null_model_vectors; i=i+2)
        {
            //clq::output("53w5", i, node_id, comm_id_neighbour);
            comm_loss += 
                internals.null_model_vectors[i][node_id]*internals.comm_loss_vectors[i+1][comm_id_neighbour]
                +internals.null_model_vectors[i+1][node_id]*internals.comm_loss_vectors[i][comm_id_neighbour];
        }
        //clq::output(comm_loss,internals.num_null_model_vectors,"loss\n");
        // gain resulting from adding to community
        double w_node_to_comm = internals.node_weight_to_communities[comm_id_neighbour];
        return markov_time * w_node_to_comm * 2 - comm_loss;
    }
};

}
