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

// Define a test fixture class template that creates the proper
// solver and auxiliar containers for the initial point, the solution and
// the bounds
template<class T>
class numerical_gradient_test : public problem_fixture<T> {
protected:
    numerical_gradient_test() = default;

    ~numerical_gradient_test() = default;
};


using testing::Types;
typedef Types<std::vector<double>, arma::vec, arma::colvec, Eigen::RowVectorXd,
        Eigen::VectorXd, std::array<double, 2> > Implementations;
TYPED_TEST_CASE(numerical_gradient_test, Implementations);

TYPED_TEST(numerical_gradient_test, rosenbrock) {
    int n = 2;
    TypeParam x, gr;
    this->set_up(std::shared_ptr<problem<TypeParam> >(new rosenbrock_function<TypeParam>(n)));
    this->random_gradient(x, gr);
    TypeParam ngr = l_bfgs_b_utils::numerical_gradient(*(this->mPb), x);
    EXPECT_RELATIVE_NEAR_VECTORS(gr, ngr);
}

TYPED_TEST(numerical_gradient_test, beale) {
    TypeParam x, gr;
    this->set_up(std::shared_ptr<problem<TypeParam> >(new beale_function<TypeParam>()));
    this->random_gradient(x, gr);
    TypeParam ngr = l_bfgs_b_utils::numerical_gradient(*(this->mPb), x);
    EXPECT_RELATIVE_NEAR_VECTORS(gr, ngr);
}


TYPED_TEST(numerical_gradient_test, booth) {
    TypeParam x, gr;
    this->set_up(std::shared_ptr<problem<TypeParam> >(new booth_function<TypeParam>()));
    this->random_gradient(x, gr);
    TypeParam ngr = l_bfgs_b_utils::numerical_gradient(*(this->mPb), x);
    EXPECT_RELATIVE_NEAR_VECTORS(gr, ngr);
}

TYPED_TEST(numerical_gradient_test, matyas) {
    TypeParam x, gr;
    this->set_up(std::shared_ptr<problem<TypeParam> >(new matyas_function<TypeParam>()));
    this->random_gradient(x, gr);
    TypeParam ngr = l_bfgs_b_utils::numerical_gradient(*(this->mPb), x);
    EXPECT_RELATIVE_NEAR_VECTORS(gr, ngr);
}

TYPED_TEST(numerical_gradient_test, goldstein) {
    TypeParam x, gr;
    this->set_up(std::shared_ptr<problem<TypeParam> >(new goldstein_price_function<TypeParam>()));
    this->random_gradient(x, gr);
    TypeParam ngr = l_bfgs_b_utils::numerical_gradient(*(this->mPb), x);
    EXPECT_RELATIVE_NEAR_VECTORS(gr, ngr);
}


