# include <vector>
# include <cmath>
# include "../cliques/math.h"
# include "../cliques/io.h"

// function to compute the matrix exponential given an input adjacency list
// comment out print statements as necessary
int main(int argc, char *argv []) {

    // create adjacency matrix in vectorised form
    std::vector<double> A = clq::read_edgelist_weighted(argv[1]);
    int size = std::sqrt(A.size());
    std::cout << "MATRIX SIZE: " << size << std::endl;

    std::cout << "INPUT MATRIX" << std::endl;
    clq::print_matrix(A,size);

    double input_time;
    if (argc < 3)
    {
        input_time = 1.0;
    } else {
        input_time = std::atof(argv[2]);
    }
    std::cout << "TIME: "<< input_time << std::endl;

    std::vector<double> result = clq::exp(A,input_time,size);
    clq::print_matrix(result,size);

    return 0;
}
