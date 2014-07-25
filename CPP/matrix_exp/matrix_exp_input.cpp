# include <cstdlib>
# include <iostream>
# include <cmath>
# include <complex>

#include "io.h"


using namespace std;

# include "matrix_exponential.hpp"
# include "c8lib.hpp"
# include "r8lib.hpp"




int main(int argc, char *argv []){

    
    // create adjacency matrix in vectorised form
    std::vector<double> A = clq::read_edgelist_weighted(argv[1]);
    int size = std::sqrt(A.size()); 
    std::cout << size << std::endl;
    
    double* Adj;
    Adj = &A[0];
    //r8mat_print(size, size, Adj, "INPUT");

    double* result = r8mat_expm1(size, Adj);
    r8mat_print(size,size,result, "RESULT");
    
    return 0;
}
