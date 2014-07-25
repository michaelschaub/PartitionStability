# include <cstdlib>
# include <iostream>
# include <cmath>
# include <complex>

using namespace std;

# include "matrix_exponential.hpp"
# include "c8lib.hpp"
# include "r8lib.hpp"

int main(){

    // Note that matrix has to be stored column first order, i.e. A(i,j) becomes A[i+N*j]
    double A[2*2]= {1.0, 1.0, 0.0, 2.0 };
    r8mat_print(2,2,A,"INPUT");

    double* result = r8mat_expm1(2, A);
    r8mat_print(2,2,result, "RESULT");
    return 0;
}
