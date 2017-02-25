#include <iostream>
#include <cmath>
#include <array>
#include "../include/problem.h"
#include "../include/l_bfgs_b.h"

using namespace std;

template<class T>
class quadratic_problem : public problem<T> {
public:
    // here, 2 is the problem input dimension
    quadratic_problem() : problem<T>(2){}

    double operator()(const T& x) {
        return std::pow(x[0] - 0.5, 2) +  std::pow(x[1] - 1, 2);
    }

    void gradient(const T& x, T& gr) {
        gr[0] = 2 * (x[0] - 0.5);
        gr[1] = 2 * (x[1] - 1);
    }
};

int main() {
    typedef std::array<double,2> my_vector;
    quadratic_problem<my_vector> qp;
    qp.set_lower_bound({-2, -2});
    qp.set_upper_bound({2, 2});

    my_vector initPoint = {2, 3};
    l_bfgs_b<my_vector> solver;
    solver.optimize(qp, initPoint);
    // Minimum should be at (0.5, 1)
    std::cout << "MINIMUM LOCATED AT : (" << initPoint[0] << ", " << initPoint[1] << ")" << std::endl;

    return 0;
}

