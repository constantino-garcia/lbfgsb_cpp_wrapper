/*
 * Copyright Constantino Antonio Garcia 2017
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "gtest/gtest.h"
#include "test_functions.h"
#include <armadillo>
#include <Eigen/Dense>
#include <initializer_list>
#include <array>
#include <vector>
#include "test_utils.h"


template <class T>
class problem_test_base : public testing::Test {
protected:
    problem_test_base(){}

    ~problem_test_base(){}

    simple_quadratic_problem<T> get_problem(int inputDimension) {
        return simple_quadratic_problem<T>(inputDimension);
    }

    simple_quadratic_problem<T> get_problem(int inputDimension, const T& lb, const T& ub) {
        return simple_quadratic_problem<T>(inputDimension, lb, ub);
    }

    // Need this method due to that Eigen does not accept init-list (to the best of my knowledge)
    T fill_proper_container(const std::initializer_list<double>& initList) {
        int n = initList.size();
        T x(n);
        auto it = initList.begin();
        for (int i = 0; i < n; i++, it++) {
            x[i] = *it;
        }
        return x;
    }
};


template<class U, std::size_t N>
class problem_test_base< std::array<U,N> > : public testing::Test {
protected:
    problem_test_base() {}

    ~problem_test_base() {}

    simple_quadratic_problem< std::array<U, N> > get_problem(int inputDimension) {
        return simple_quadratic_problem< std::array<U, N> >(inputDimension);
    }

    simple_quadratic_problem< std::array<U, N> > get_problem(int inputDimension,
                                                             const std::array<U, N>& lb,
                                                             const std::array<U, N>& ub) {
        return simple_quadratic_problem< std::array<U,N> >(inputDimension, lb, ub);
    }
    std::array<U, N> fill_proper_container(const std::initializer_list<double>& initList) {
        assert(N == initList.size());
        std::array<U, N> x;
        std::copy(initList.begin(), initList.end(), std::begin(x));
        return x;
    }

};

template <class T>
class problem_test : public problem_test_base<T> {};
template <class T>
class dyn_problem_test : public problem_test_base<T> {};
// don't test Eigen vectors to make use of initializer_lists

using testing::Types;
typedef Types<std::vector<double>, arma::vec, std::array<double, 3> > Implementations;
TYPED_TEST_CASE(problem_test, Implementations);

TYPED_TEST(problem_test, constructor_invalid_dim) {
    EXPECT_THROW(this->get_problem(0), std::invalid_argument);
    EXPECT_THROW(this->get_problem(-1), std::invalid_argument);
}


TYPED_TEST(problem_test, constructor_valid_dim) {
    simple_quadratic_problem<TypeParam> pb= this->get_problem(3);
    EXPECT_EQ(3, pb.getInputDimension());
}


TYPED_TEST(problem_test, constructor_valid_bounds) {
    TypeParam lb = {1, 2, 3};
    TypeParam ub = {4, 5, 6};
    simple_quadratic_problem<TypeParam> pb(3, lb, ub);

    EXPECT_EQ_VECTORS(lb, pb.get_lower_bound());
    EXPECT_EQ_VECTORS(ub, pb.get_upper_bound());
}


TYPED_TEST(problem_test, constructor_invalid_bounds) {
    TypeParam lb1 = {5, 2, 3};
    TypeParam lb2 = {1, 6, 3};
    TypeParam lb3 = {1, 2, 7};
    TypeParam ub = {4, 5, 6};
    EXPECT_THROW(this->get_problem(3, lb1, ub), std::invalid_argument);
    EXPECT_THROW(this->get_problem(3, lb2, ub), std::invalid_argument);
    EXPECT_THROW(this->get_problem(3, lb3, ub), std::invalid_argument);
}


TYPED_TEST(problem_test, setter_valid_bound) {
    TypeParam lb = {3, 4, 5};
    TypeParam ub = {3.5, 4.5, 5.5};
    simple_quadratic_problem<TypeParam> pb = this->get_problem(3, {1, 2, 3}, {4, 5, 6});
    pb.set_lower_bound(lb);
    pb.set_upper_bound(ub);

    EXPECT_EQ_VECTORS(lb, pb.get_lower_bound());
    EXPECT_EQ_VECTORS(ub, pb.get_upper_bound());
}




TYPED_TEST(problem_test, setter_invalid_bounds) {
    simple_quadratic_problem<TypeParam> pb = this->get_problem(3, {1, 2, 3}, {4, 5, 6});
    TypeParam lb1 = {5, 2, 3};
    TypeParam lb2 = {1, 6, 3};
    TypeParam lb3 = {1, 2, 7};
    TypeParam ub1 = {0, 5, 6};
    TypeParam ub2 = {1, 1, 6};
    TypeParam ub3 = {3, 5, 2};
    EXPECT_THROW(pb.set_lower_bound(lb1), std::invalid_argument);
    EXPECT_THROW(pb.set_lower_bound(lb2), std::invalid_argument);
    EXPECT_THROW(pb.set_lower_bound(lb3), std::invalid_argument);
    EXPECT_THROW(pb.set_upper_bound(ub1), std::invalid_argument);
    EXPECT_THROW(pb.set_upper_bound(ub2), std::invalid_argument);
    EXPECT_THROW(pb.set_upper_bound(ub3), std::invalid_argument);
}


// these tests are not applicable for the std::array-based problems, since the
// consistency of the arrays are check in compilation time
typedef Types<std::vector<double>, arma::vec > dynImplementations;
TYPED_TEST_CASE(dyn_problem_test, dynImplementations);


TYPED_TEST(dyn_problem_test, constructor_invalid_sized_bounds) {
    TypeParam  b = {4, 5, 6};
    EXPECT_THROW(this->get_problem(3, {1, 2}, b), std::invalid_argument);
    EXPECT_THROW(this->get_problem(3, b, {1, 2}), std::invalid_argument);
    EXPECT_THROW(this->get_problem(3, {1, 2, 3, 4}, b), std::invalid_argument);
    EXPECT_THROW(this->get_problem(3, b, {1, 2, 3, 4}), std::invalid_argument);
}

TYPED_TEST(dyn_problem_test, setter_invalid_sized_bounds) {
    TypeParam b = {1, 2};
    TypeParam b2 = {1, 2, 3, 4};
    simple_quadratic_problem<TypeParam> pb = this->get_problem(3, {1, 2, 3}, {4, 5, 6});
    EXPECT_THROW(pb.set_lower_bound(b), std::invalid_argument);
    EXPECT_THROW(pb.set_lower_bound(b2), std::invalid_argument);
}
