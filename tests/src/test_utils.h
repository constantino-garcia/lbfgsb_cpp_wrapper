/*
 * Copyright Constantino Antonio Garcia 2017
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef LBFGSB_CPP_WRAPPER_TEST_UTILS_H
#define LBFGSB_CPP_WRAPPER_TEST_UTILS_H

template <class T>
inline void EXPECT_NEAR_VECTORS(const T& x, const T& y, double epsilon = 1e-5) {
    ASSERT_EQ(x.size(),y.size());
    for (int i = 0; i < x.size(); ++i) {
        ASSERT_NEAR(x[i], y[i], epsilon);
    }
}

template <class T>
inline void EXPECT_EQ_VECTORS(const T& x, const T& y) {
    ASSERT_EQ(x.size(),y.size());
    for (int i = 0; i < x.size(); ++i) {
        ASSERT_EQ(x[i], y[i]);
    }
}

#endif //LBFGSB_CPP_WRAPPER_TEST_UTILS_H
