/*
 * Copyright Constantino Antonio Garcia 2017
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef LBFGSB_CPP_UTILS_H
#define LBFGSB_CPP_UTILS_H

#include <initializer_list>
#include <cassert>

namespace l_bfgs_b_utils {
    template<class T, typename F>
    T numerical_gradient(F &functor, const T &x, const T &lowerBound, const T &upperBound,
                         double gridSpacing = 1e-6) {
        // check consistency of the dimensions
        int inputDimension = x.size();
        if (inputDimension != lowerBound.size() || inputDimension != upperBound.size()) {
            throw std::invalid_argument("The size of x does not match the bound's dimensions");
        }
        // check x is within bounds
        for (int i = 0; i < inputDimension; i++) {
            if (x[i] > upperBound[i] | x[i] < lowerBound[i]) {
                throw std::runtime_error("x is not contained within [lowerBound, upperBound]");
            }
        }
        T gr(x);
        T workX(x);
        for (int i = 0; i < inputDimension; i++) {
            double effectiveGridOver =
                    ((x[i] + gridSpacing) > upperBound[i]) ? upperBound[i] - x[i] : gridSpacing;
            workX[i] = x[i] + effectiveGridOver;
            double valueOver = functor(workX);
            double effectiveGridBelow =
                    ((x[i] - gridSpacing) < lowerBound[i]) ? x[i] - lowerBound[i] : gridSpacing;
            workX[i] = x[i] - effectiveGridBelow;
            gr[i] = (valueOver - functor(workX)) / (effectiveGridOver + effectiveGridBelow);
            // restore original value
            workX[i] = x[i];
        }
        return gr;
    }

    template<class T, typename F>
    T numerical_gradient(F &functor, const T &x, double gridSpacing = 1e-6) {
        T lowerBound(x), upperBound(x);
        int inputDimension = x.size();
        for (int i = 0; i < inputDimension; i++) {
            lowerBound[i] = -std::numeric_limits<double>::infinity();
            upperBound[i] = std::numeric_limits<double>::infinity();
        }
        return numerical_gradient(functor, x, lowerBound, upperBound, gridSpacing);
    }

    template<class T>
    void fill_container(T& x, const std::initializer_list<double> &initList) {
        int n = initList.size();
        x.resize(n);
        auto it = initList.begin();
        for (int i = 0; i < n; i++, it++) {
            x[i] = *it;
        }
    }

    template<std::size_t N>
    void fill_container(std::array<double, N>& x, const std::initializer_list<double>& initList) {
        assert(N == initList.size());
        std::copy(initList.begin(), initList.end(), std::begin(x));
    }
}
#endif //LBFGSB_CPP_UTILS_H