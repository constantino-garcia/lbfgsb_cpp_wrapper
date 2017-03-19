
#ifndef LBFGSB_CPP_PROBLEM_FIXTURE_H_H
#define LBFGSB_CPP_PROBLEM_FIXTURE_H_H

#include <lbfgsb_cpp/problem.h>
#include <memory>
#include "random_vector_generator.h"

template<class T>
class problem_fixture : public testing::Test {
public:
    problem_fixture() = default;

    ~problem_fixture() = default;

    void set_up(std::shared_ptr<problem<T> > ptr) {
        mPb = ptr;
        // get range for generating random vectors from the problem's bounds
        T lb = mPb->get_lower_bound();
        T ub = mPb->get_upper_bound();
        double minValue = lb[0];
        double maxValue = ub[0];
        // get most restrictive values from the upper and lower bounds
        for (int i = 1; i < mPb->get_input_dimension(); i++) {
            if (lb[i] > minValue) {
                minValue = lb[i];
            }
            if (ub[i] < maxValue) {
                maxValue = ub[i];
            }
        }
        // The random generator does not handle properly infinite values:
        // eliminate them
        if (std::isinf(maxValue) && std::isinf(minValue)) {
            minValue = -100;
            maxValue = 100;
        } else if (std::isinf(maxValue)) {
            maxValue = minValue + 100;
        } else if (std::isinf(minValue)) {
            minValue = maxValue - 100;
        }
        mRvg = random_vector_generator<T>(mPb->get_input_dimension(),
                                          minValue, maxValue);

    }

    void random_point(T &x) {
        x = mRvg();
    }

    // calculate gradient at a random point
    void random_gradient(T &x, T &gr) {
        x = mRvg();
        // ensure the proper dimensions of gr without relying on the container constructor
        gr = mRvg();
        mPb->gradient(x, gr);
    }

protected:
    std::shared_ptr<problem<T> > mPb;
    T pbMinimum;
    random_vector_generator<T> mRvg;
};

#endif //LBFGSB_CPP_PROBLEM_FIXTURE_H_H
