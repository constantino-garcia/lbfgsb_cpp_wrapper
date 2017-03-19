/*
 * Copyright Constantino Antonio Garcia 2017
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "gtest/gtest.h"
#include "test_functions.h"
#include "test_utils.h"
#include "problem_fixture.h"
#include <lbfgsb_cpp/l_bfgs_b.h>
#include <lbfgsb_cpp/utils.h>
#include <vector>
#include <armadillo>
#include <Eigen/Dense>
#include <array>

template<class T>
class l_bfgs_b_num_gradient_test : public problem_fixture<T> {
public:
    l_bfgs_b_num_gradient_test() {
        mSolver.set_machine_precision_factor(10);
        mSolver.set_projected_gradient_tolerance(0);
        mSolver.set_max_iterations(1000);
    }

    ~l_bfgs_b_num_gradient_test() = default;

    void test_optimization(const std::initializer_list<double> &solution) {
        T x, sol;
        l_bfgs_b_utils::fill_container(sol, solution);
        for (int i = 0; i < mNoTests; i++) {
            this->random_point(x);
            this->mSolver.optimize(*(this->mPb), x);
            EXPECT_NEAR_VECTORS(sol, x, mTolerance);
        }
    }

protected:
    l_bfgs_b<T> mSolver;
    int mNoTests = 100;
    double mTolerance = 1e-3;
};


using testing::Types;
typedef Types<std::vector<double>, arma::vec, arma::colvec, Eigen::RowVectorXd,
        Eigen::VectorXd, std::array<double, 2> > Implementations;
TYPED_TEST_CASE(l_bfgs_b_num_gradient_test, Implementations);

TYPED_TEST(l_bfgs_b_num_gradient_test, rosenbrock) {
    int n = 2;
    std::shared_ptr<problem<TypeParam> > ptr(new rosenbrock_function<TypeParam>(n));
    ptr->set_lower_bound({-10, -10});
    ptr->set_upper_bound({10, 10});
    this->set_up(ptr);
    this->test_optimization({1, 1});
}

TYPED_TEST(l_bfgs_b_num_gradient_test, rosenbrock_numerical_gradient) {
    int n = 2;
    std::shared_ptr<problem<TypeParam> > ptr(new rosenbrock_function_base<TypeParam>(n));
    ptr->set_lower_bound({-10, -10});
    ptr->set_upper_bound({10, 10});
    this->set_up(ptr);
    this->test_optimization({1, 1});
}

// some initial points of beale function do not converge to the global minimum
// (tested with R's L-BFGS_B implementations). Restrict the searching domain to
// avoid this issue.
TYPED_TEST(l_bfgs_b_num_gradient_test, beale) {
    std::shared_ptr<problem<TypeParam> > ptr(new beale_function<TypeParam>());
    ptr->set_lower_bound({0, -2});
    ptr->set_upper_bound({4.5, 1});
    this->set_up(ptr);
    this->test_optimization({3, 0.5});
}

TYPED_TEST(l_bfgs_b_num_gradient_test, beale_numerical_gradient) {
    std::shared_ptr<problem<TypeParam> > ptr(new beale_function_base<TypeParam>());
    ptr->set_lower_bound({0, -2});
    ptr->set_upper_bound({4.5, 1});
    this->set_up(ptr);
    this->test_optimization({3, 0.5});
}

TYPED_TEST(l_bfgs_b_num_gradient_test, booth) {
    std::shared_ptr<problem<TypeParam> > ptr(new booth_function<TypeParam>());
    ptr->set_lower_bound({-10, -10});
    ptr->set_upper_bound({10, 10});
    this->set_up(ptr);
    this->test_optimization({1, 3});
}

TYPED_TEST(l_bfgs_b_num_gradient_test, booth_numerical_gradient) {
    std::shared_ptr<problem<TypeParam> > ptr(new booth_function_base<TypeParam>());
    ptr->set_lower_bound({-10, -10});
    ptr->set_upper_bound({10, 10});
    this->set_up(ptr);
    this->test_optimization({1, 3});
}

TYPED_TEST(l_bfgs_b_num_gradient_test, matyas) {
    std::shared_ptr<problem<TypeParam> > ptr(new matyas_function<TypeParam>());
    ptr->set_lower_bound({-10, -10});
    ptr->set_upper_bound({10, 10});
    this->set_up(ptr);
    this->test_optimization({0, 0});
}

TYPED_TEST(l_bfgs_b_num_gradient_test, matyas_numerical_gradient) {
    std::shared_ptr<problem<TypeParam> > ptr(new matyas_function_base<TypeParam>());
    ptr->set_lower_bound({-10, -10});
    ptr->set_upper_bound({10, 10});
    this->set_up(ptr);
    this->test_optimization({0, 0});
}


// goldstein function has several local minima: restrict the search space to find the
// global minimum
TYPED_TEST(l_bfgs_b_num_gradient_test, goldstein) {
    std::shared_ptr<problem<TypeParam> > ptr(new goldstein_price_function<TypeParam>());
    ptr->set_lower_bound({-2, -2});
    ptr->set_upper_bound({2, -0.75});
    this->set_up(ptr);
    this->test_optimization({0, -1});
}

TYPED_TEST(l_bfgs_b_num_gradient_test, goldstein_numerical_gradient) {
    std::shared_ptr<problem<TypeParam> > ptr(new goldstein_price_function_base<TypeParam>());
    ptr->set_lower_bound({-2, -2});
    ptr->set_upper_bound({2, -0.75});
    this->set_up(ptr);
    this->test_optimization({0, -1});
}

