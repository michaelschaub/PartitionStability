#pragma once

#include "linearised_internals_generic.h"
#include "stability_gen.h"

namespace clq {
/////////////////////////////////////////////
// NORMALISED STABILIY
/////////////////////////////////////////////
// Linearised normalised stability w/o partition
//template<typename G, typename M>
//clq::LinearisedInternals gen(find_linearised_normalised_stability,
//		G &graph, M &weights) {
//	LinearisedInternals internals(graph, weights);
//	return internals;
//}
//
// Linearised normalised stability with partition given
//template<typename G, typename M, typename P>
//clq::LinearisedInternals gen(find_linearised_normalised_stability,
//		G &graph, M &weights, P &partition, P partition_unused = P(0),
//		std::vector<double> null_model_vec = std::vector<double>()) {
//	LinearisedInternals internals(graph, weights, partition);
//	return internals;
//}

/////////////////////////////////////////////
// Generalized Internals
/////////////////////////////////////////////
template<typename G, typename M>
clq::LinearisedInternalsGeneric gen(find_linearised_generic_stability,
		G &graph, M &weights, std::vector<std::vector<double>> null_model) {
	LinearisedInternalsGeneric internals(graph, weights, null_model);
	return internals;
}

template<typename G, typename M, typename P>
clq::LinearisedInternalsGeneric gen(find_linearised_generic_stability,
		G &graph, M &weights, P &partition, P partition_unused = P(0), std::vector<std::vector<double>> null_model = std::vector<std::vector<double>>()) {
	LinearisedInternalsGeneric internals(graph, weights, partition, null_model);
	return internals;
}

}
