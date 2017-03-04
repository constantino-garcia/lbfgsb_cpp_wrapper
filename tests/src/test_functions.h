/*
 * Copyright Constantino Antonio Garcia 2017
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef L_BFGS_B_CPP_WRAPPER_TEST_FUNCTIONS_H
#define L_BFGS_B_CPP_WRAPPER_TEST_FUNCTIONS_H

#include <lbfgsb_cpp/problem.h>
#include <algorithm>


// See https://en.wikipedia.org/wiki/Test_functions_for_optimization for a complete list of numerical optimization tests

template<class T>
class rosenbrock_function : public problem<T> {
public:
    rosenbrock_function(int n) : problem<T>(n) {};

    virtual ~rosenbrock_function() = default;

    double operator()(const T &x) {
        double result = 0.0;
        for (int i = 0; i < (this->mInputDimension - 1); ++i) {
            result += (100 * std::pow(x[i + 1] - std::pow(x[i], 2), 2.0) + std::pow(1 - x[i], 2.0));
        }
        return result;
    }

    void gradient(const T &x, T &gr) {
        if (gr.size() != this->mInputDimension) {
            throw std::invalid_argument("gradient size does not match input dimension");
        }
        // resize the gradient to avoid positive increments of the functions
        gr[0] = (400 * (std::pow(x[0], 3.0) - x[1] * x[0]) + 2 * (x[0] - 1));
        for (int i = 1; i < (this->mInputDimension - 1); ++i) {
            gr[i] = (200 * (x[i] - std::pow(x[i - 1], 2.0)) +
                     400 * (std::pow(x[i], 3.0) - x[i + 1] * x[i]) +
                     2 * (x[i] - 1));
        }
        gr[this->mInputDimension - 1] =
                (200 * (x[this->mInputDimension - 1] - std::pow(x[this->mInputDimension - 2], 2.0)));
    }
};


template<class T>
class beale_function : public problem<T> {
public:
    beale_function() : problem<T>(2) {};

    virtual ~beale_function() = default;

    double operator()(const T &x) {
        return std::pow(1.5 - x[0] + x[0] * x[1], 2) +
               std::pow(2.25 - x[0] + x[0] * std::pow(x[1], 2), 2) +
               std::pow(2.625 - x[0] + x[0] * std::pow(x[1], 3), 2);
    }

    void gradient(const T &x, T &gr) {
        if (gr.size() != this->mInputDimension) {
            throw std::invalid_argument("gradient size does not match input dimension");
        }
        gr[0] = (2 * (1.5 - x[0] + x[0] * x[1]) * (-1 + x[1]) +
                 2 * (2.25 - x[0] + x[0] * std::pow(x[1], 2)) * (-1 + std::pow(x[1], 2)) +
                 2 * (2.625 - x[0] + x[0] * std::pow(x[1], 3)) * (-1 + std::pow(x[1], 3)));
        gr[1] = (2 * (1.5 - x[0] + x[0] * x[1]) * x[0] +
                 2 * (2.25 - x[0] + x[0] * std::pow(x[1], 2)) * 2 * x[0] * x[1] +
                 2 * (2.625 - x[0] + x[0] * std::pow(x[1], 3)) * 3 * x[0] * std::pow(x[1], 2));
    }
};


template<class T>
class goldstein_price_function : public problem<T> {
public:
    goldstein_price_function() : problem<T>(2) {};

    virtual ~goldstein_price_function() = default;

    double operator()(const T &x) {
        return (1 + std::pow(x[0] + x[1] + 1, 2) *
                    (19 - 14 * x[0] + 3 * x[0] * x[0] - 14 * x[1] + 6 * x[0] * x[1] + 3 * x[1] * x[1])) *
               (30 + std::pow(2 * x[0] - 3 * x[1], 2) *
                     (18 - 32 * x[0] + 12 * x[0] * x[0] + 48 * x[1] - 36 * x[0] * x[1] + 27 * x[1] * x[1]));
    }

    void gradient(const T &xx, T &gr) {
        if (gr.size() != this->mInputDimension) {
            throw std::invalid_argument("gradient size does not match input dimension");
        }
        double x = xx[0];
        double y = xx[1];
        // as output by wolfram alpha
        gr[0] = 24 * (8 * x * x * x - 4 * x * x * (9 * y + 4) + 6 * x * (9 * y * y + 8 * y + 1) -
                      9 * y * (3 * y * y + 4 * y + 1)) * ((3 * x * x + 2 * x * (3 * y - 7) + 3 * y * y - 14 * y + 19) *
                                                          std::pow(x + y + 1, 2) + 1) +
                12 * (x * x * x + x * x * (3 * y - 2) + x * (3 * y * y - 4 * y - 1) + y * y * y -
                      2 * y * y - y + 2) * ((12 * x * x - 4 * x * (9 * y + 8) + 3 * (9 * y * y + 16 * y + 6)) *
                                            std::pow(2 * x - 3 * y, 2) + 30);
        gr[1] = 12 * (x * x * x + x * x * (3 * y - 2) + x * (3 * y * y - 4 * y - 1) + y * y * y - 2 * y * y - y + 2) *
                ((12 * x * x - 4 * x * (9 * y + 8) + 3 * (9 * y * y + 16 * y + 6)) * std::pow(2 * x - 3 * y, 2) + 30) -
                36 * (8 * x * x * x - 4 * x * x * (9 * y + 4) + 6 * x * (9 * y * y + 8 * y + 1) - 9 * y * (3 * y * y +
                                                                                                           4 * y + 1)) *
                ((3 * x * x + 2 * x * (3 * y - 7) + 3 * y * y - 14 * y + 19) * std::pow(x + y + 1, 2) + 1);
    }
};

template<class T>
class booth_function : public problem<T> {
public:
    booth_function() : problem<T>(2) {};

    virtual ~booth_function() = default;

    double operator()(const T &x) {
        return std::pow(x[0] + 2 * x[1] - 7, 2) + std::pow(2 * x[0] + x[1] - 5, 2);
    }

    void gradient(const T &x, T &gr) {
        if (gr.size() != this->mInputDimension) {
            throw std::invalid_argument("gradient size does not match input dimension");
        }
        gr[0] = (2 * (x[0] + 2 * x[1] - 7) + 4 * (2 * x[0] + x[1] - 5));
        gr[1] = (4 * (x[0] + 2 * x[1] - 7) + 2 * (2 * x[0] + x[1] - 5));
    }
};


template<class T>
class matyas_function : public problem<T> {
public:
    matyas_function() : problem<T>(2) {};

    virtual ~matyas_function() = default;

    double operator()(const T &x) {
        return 0.26 * (x[0] * x[0] + x[1] * x[1]) - 0.48 * x[0] * x[1];
    }

    void gradient(const T &x, T &gr) {
        if (gr.size() != this->mInputDimension) {
            throw std::invalid_argument("gradient size does not match input dimension");
        }
        gr[0] = 0.26 * 2 * x[0] - 0.48 * x[1];
        gr[1] = 0.26 * 2 * x[1] - 0.48 * x[0];
    }
};


template<class T>
class simple_quadratic_problem : public problem<T> {
public:
    simple_quadratic_problem(int inputDimension) : problem<T>(inputDimension) {}

    simple_quadratic_problem(int inputDimension, const T &lb, const T &ub)
            : problem<T>(inputDimension, lb, ub) {}

    double operator()(const T &x) {
        return std::accumulate(std::begin(x), std::end(x), 0,
                               [](double acc, double value) {
                                   return acc + value * value;
                               });
    }

    void gradient(const T &x, T &gr) {
        std::transform(x.begin(), x.end(), gr.begin(),
                       [](double value) { return 2 * value; });
    }

};

#endif //L_BFGS_B_CPP_WRAPPER_TEST_FUNCTIONS_H

