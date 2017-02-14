#ifndef L_BFGS_B_CPP_WRAPPER_TEST_FUNCTIONS_H
#define L_BFGS_B_CPP_WRAPPER_TEST_FUNCTIONS_H

#include "../../include/problem.h"

template <class T>
class rosenbrock_function : public problem<T> {
public:
    rosenbrock_function():  problem<T>(2) {};

    virtual ~rosenbrock_function() = default;

    double operator()(const T& x) {
        double x1 = x[0];
        double x2 = x[1];

        return (100 * std::pow(x2 - std::pow(x1, 2), 2.0) + std::pow(1 - x1, 2.0));
    }

    void gradient(const T& x, T& gr) {
        double x1 = x[0];
        double x2 = x[1];

        double epsilon = 1e-3;
        gr.set_size(2);
        gr(0) = epsilon * (-2 * (1 - x1) + 400 * (std::pow(x1, 3.0) - x2 * x1));
        gr(1) = epsilon * (200 * (x2 - std::pow(x1, 2.0)));

    }
};

#endif //L_BFGS_B_CPP_WRAPPER_TEST_FUNCTIONS_H

