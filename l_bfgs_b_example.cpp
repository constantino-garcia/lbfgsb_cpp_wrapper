#include <iostream>
#include <armadillo>
#include "tests/src/test_functions.h"
#include "l_bfgs_b.h"

using namespace std;

int main() {
    rosenbrock_function<arma::vec> rf;
    rf.set_lower_bound({-10,-10});
    rf.set_upper_bound({10,10});

    l_bfgs_b<arma::vec> solver(5);
    solver.setVerboseLevel(1);
    arma::vec x0 = {-1, 2};
    solver.optimize(rf, x0);

    std::cout << "***************************" << std::endl;
    std::cout << "MINIMUM LOCATED AT : (" << x0[0] << ", " << x0[1] << ")"<< std::endl;

    return 0;
}

