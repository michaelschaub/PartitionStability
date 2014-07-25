# include <vector>
# include "../cliques/math.h"
# include "../cliques/io.h"


int main(){

    // Note that matrix has to be stored column first order, i.e. A(i,j) becomes A[i+N*j]
    double a[2*2]= {1.0, 1.0, 0.0, 2.0 };
    std::vector<double> A (a,a+sizeof(a)/sizeof(double) );
    clq::print_matrix(A,2);
    
    double time = 1.0;

    std::vector<double> result = clq::exp(A,time,2);
    clq::print_matrix(result,2);
    return 0;
}
