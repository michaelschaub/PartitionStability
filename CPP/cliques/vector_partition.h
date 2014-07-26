#pragma once

#include <vector>
#include <set>
#include <map>
#include <iostream>

namespace clq {
/**
 @brief  A partition where the position in the vector denotes the node_id
 and its value determines its set id.

 Good when you need constant time assignment of a node to a set
 and removal of a node from a set
 */
class VectorPartition {
private:
	int num_nodes;
	std::vector<int> partition_vector;
	bool is_normalised;

public:
	//#################### CONSTRUCTORS ####################

	// construct empty partition
	explicit VectorPartition(int num_nodes) :
		num_nodes(num_nodes),
				partition_vector(std::vector<int>(num_nodes, -1)),
				is_normalised(false) {
	}

	// construct partition with initial_set
	explicit VectorPartition(int num_nodes, int initial_set) :
		num_nodes(num_nodes), partition_vector(std::vector<int>(num_nodes,
				initial_set)), is_normalised(false) {
	}

	// construct partition from vector
	explicit VectorPartition(std::vector<int> partition) :
		num_nodes(partition.size()), partition_vector(partition),
				is_normalised(false) {
	}
	// what is that for a horribly over complicated syntax?? --M
	void initialise_as_singletons() {
		for (std::vector<int>::iterator itr = partition_vector.begin(); itr
				!= partition_vector.end(); ++itr) {
			*itr = itr - partition_vector.begin();
		}
		is_normalised = true;
	}

	void initialise_as_global() {
		partition_vector = std::vector<int>(partition_vector.size(), 0);
	}

	//#################### PUBLIC METHODS ####################
public:
	int find_set(int node_id) const {
		return partition_vector[node_id];
	}

	void isolate_node(int node_id) {
		partition_vector[node_id] = -1;
		is_normalised = false;
	}

	void add_node_to_set(int node_id, int set_id) {
		partition_vector[node_id] = set_id;
		is_normalised = false;
	}

	// Modifies ids such that they are contiguous, start at 0
	// e.g. 2,1,4,2 -> 0,1,2,0
	std::map<int, int> normalise_ids() {
		// A mapping from new set ids to old set ids
		// May be useful in some algorithms
		std::map<int, int> set_new_to_old;
		if (!is_normalised) {
			int start_num = 0;
			std::map<int, int> set_old_to_new;

			// For every element, make a map
			for (std::vector<int>::iterator itr = partition_vector.begin(); itr
					!= partition_vector.end(); ++itr) {

				// Map the old set_id to a new one (starting at 0)
				// Unless the mapping has already been done (then look it up)
				std::map<int, int>::iterator old_set =
						set_old_to_new.find(*itr);
				if (old_set == set_old_to_new.end()) {
					set_old_to_new[*itr] = start_num;
					set_new_to_old[start_num] = *itr;
					*itr = start_num;
					start_num++;
				} else {
					*itr = old_set->second;
				}
			}
			is_normalised = true;
			return set_new_to_old;
		} else {
			// Still need to reconstruct new to old mapping (since we don't store it)
			for (std::vector<int>::iterator itr = partition_vector.begin(); itr
					!= partition_vector.end(); ++itr) {
				set_new_to_old[*itr] = *itr;
			}
			return set_new_to_old;
		}
	}

    void normalise_id_vector(){
        // create tempory vector to store mapping and a counter variable for the new ids assigned so far
        std::vector<int> old_to_new(num_nodes,-1);
        int new_id_assignment = 0;

        if (!is_normalised) {
	        for(int i=0; i<num_nodes; i++){
    	        // current, i.e. "old" id of node
    	        int old_id = partition_vector[i];

    	        // Unassigned node (e.g. from isolate node)
    	        if (old_id == -1) {
    	        	partition_vector[i] = new_id_assignment;
    	        	new_id_assignment++;
    	        }

	            // if the old id has not been mapped to a new one so far, 
	            // map it to a new id and adjust partition vector
	            else if(old_to_new[old_id] == -1){
	                old_to_new[old_id] = new_id_assignment;
	                partition_vector[i] = new_id_assignment;
	                new_id_assignment++;
	            }
	            // otherwise apply map
	            else{
	                partition_vector[i] = old_to_new[old_id];
	            }
	        }
	        is_normalised = true;
    	}
    }

	int element_count() const {
		return partition_vector.size();
	}

	int set_count() {
		std::set<int> seen_nodes;
		for (std::vector<int>::iterator itr = partition_vector.begin(); itr
				!= partition_vector.end(); ++itr) {
			if (*itr != -1) {
				seen_nodes.insert(*itr);
			}
		}
		return seen_nodes.size();
	}

	std::vector<int> get_nodes_from_set(int set_id) {
		std::vector<int> nodes_in_set;
		for (int i = 0; i < num_nodes; ++i) {
			if (partition_vector[i] == set_id) {
				nodes_in_set.push_back(i);
			}
		}
		return nodes_in_set;
	}

	std::vector<int> return_partition_vector() {
		return partition_vector;
	}

	//#################### Operators ####################
	bool operator==(const VectorPartition& other) const {
		if (this->is_normalised && other.is_normalised) {
			return (other.partition_vector == this->partition_vector);
		} else {
			VectorPartition a(*this);
			VectorPartition b(other);

			a.normalise_ids();
			b.normalise_ids();
			return a.partition_vector == b.partition_vector;
		}
	}

	bool operator!=(const VectorPartition& other) const {
		return !(partition_vector == other.partition_vector);
	}
};

}
