/*
 * Copyright Constantino Antonio Garcia 2017
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef LBFGSB_CPP_TEST_UTILS_H
#define LBFGSB_CPP_TEST_UTILS_H

#include "gtest/gtest.h"
#include <cassert>

template<class T>
inline void EXPECT_NEAR_VECTORS(const T &x, const T &y, double epsilon = 1e-5) {
    ASSERT_EQ(x.size(), y.size());
    for (int i = 0; i < x.size(); ++i) {
        ASSERT_NEAR(x[i], y[i], epsilon);
    }
}

// use relative error
template<class T>
inline void EXPECT_RELATIVE_NEAR_VECTORS(const T &x, const T &y, double epsilon = 1e-3) {
    ASSERT_EQ(x.size(), y.size());
    for (int i = 0; i < x.size(); ++i) {
        ASSERT_NEAR((x[i] - y[i]) / x[i], 0, epsilon);
    }
}

template<class T>
inline void EXPECT_EQ_VECTORS(const T &x, const T &y) {
    ASSERT_EQ(x.size(), y.size());
    for (int i = 0; i < x.size(); ++i) {
        ASSERT_EQ(x[i], y[i]);
    }
}


// Need this function due to that Eigen does not accept init-list (to the best of my knowledge)
// The size of x should match the size of initializer_list
template<class T>
void fill_proper_container(const std::initializer_list<double> &initList, T &x) {
    int n = initList.size();
    assert(n == x.size());
    auto it = initList.begin();
    for (int i = 0; i < n; i++, it++) {
        x[i] = *it;
    }
}

template<class U, std::size_t N>
void fill_proper_container(const std::initializer_list<double> &initList, std::array<U, N> &x) {
    assert(N == initList.size());
    std::copy(initList.begin(), initList.end(), std::begin(x));
}

#endif //LBFGSB_CPP_TEST_UTILS_H
