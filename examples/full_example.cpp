/*
 * Copyright Constantino Antonio Garcia 2017
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <iostream>
#include <armadillo>
#include <Eigen/Dense>
#include "../tests/src/test_functions.h"
#include <lbfgsb_cpp/problem.h>
#include <lbfgsb_cpp/l_bfgs_b.h>


using namespace std;

int main() {
    std::cout << "*** MINIMIZING BOOTH (with std::vector) ***" << std::endl;
    booth_function<std::vector<double> > bl;
    bl.set_lower_bound({-2, -2});
    bl.set_upper_bound({5, 4});
    l_bfgs_b<std::vector<double> > stdSolver(5);
    stdSolver.set_verbose_level(-1);
    std::vector<double> stdX = {-1, 2};
    stdSolver.optimize(bl, stdX);
    std::cout << "MINIMUM LOCATED AT : (" << stdX[0] << ", " << stdX[1] << ")" << std::endl << std::endl;

    std::cout << "*** MINIMIZING BEALE (with std::array) ***" << std::endl;
    beale_function<std::array<double, 2> > beale;
    beale.set_lower_bound({-4.5, -4.5});
    beale.set_upper_bound({4.5, 4.5});
    l_bfgs_b<std::array<double, 2> > arraySolver(5);
    arraySolver.set_verbose_level(-1);
    // The Beale function has a explosive gradient near the corners... scale the gradient
    arraySolver.set_gradient_scaling_factor(1e-2);
    std::array<double, 2> arrayX = {-1, 0};
    arraySolver.optimize(beale, arrayX);
    std::cout << "MINIMUM LOCATED AT : (" << arrayX[0] << ", " << arrayX[1] << ")" << std::endl << std::endl;

    std::cout << "*** MINIMIZING GOLDSTEIN (with armadillo) ***" << std::endl;
    goldstein_price_function<arma::vec> gold;
    gold.set_lower_bound({-2, -2});
    gold.set_upper_bound({2, 2});
    l_bfgs_b<arma::vec> armaSolver(5);
    armaSolver.set_verbose_level(-1);
    arma::vec armaX = {-1, 1.5};
    armaSolver.optimize(gold, armaX);
    std::cout << "MINIMUM LOCATED AT : (" << armaX[0] << ", " << armaX[1] << ")" << std::endl << std::endl;

    std::cout << "*** MINIMIZING MATYAS (with Eigen3) ***" << std::endl;
    matyas_function<Eigen::VectorXd> matyas;
    Eigen::VectorXd eigenX(2), lb(2), ub(2);
    lb << -10, -10;
    ub << 10, 10;
    matyas.set_lower_bound(lb);
    matyas.set_upper_bound(ub);
    l_bfgs_b<Eigen::VectorXd> eigenSolver(5);
    eigenSolver.set_verbose_level(-1);
    eigenX << 8, -8;
    eigenSolver.optimize(matyas, eigenX);
    std::cout << "MINIMUM LOCATED AT : (" << eigenX[0] << ", " << eigenX[1] << ")" << std::endl << std::endl;

    // It is possible to modify the L-BFGS-B parameters
    std::cout << "*** MINIMIZING ROSENBROCK ***" << std::endl;
    armaSolver.set_verbose_level(0);
    armaSolver.set_machine_precision_factor(1e2);
    rosenbrock_function<arma::vec> rf(3);
    rf.set_lower_bound({-10, -10, -10});
    rf.set_upper_bound({10, 10, 10});
    armaX = {-1, 2, 2};
    armaSolver.optimize(rf, armaX);
    std::cout << "***************************" << std::endl;
    std::cout << "MINIMUM LOCATED AT : (" << armaX[0] << ", " << armaX[1] << ", " << armaX[2] << ")" << std::endl;

    return 0;
}

