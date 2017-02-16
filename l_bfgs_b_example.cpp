#include <iostream>
#include <armadillo>
#include "tests/src/test_functions.h"
#include "l_bfgs_b.h"
#include <Eigen/Dense>

using namespace std;

int main() {
    l_bfgs_b<arma::vec> armaSolver(5);
    armaSolver.setVerboseLevel(-1);
    arma::vec armaX;

    std::vector<double> stdX;
    l_bfgs_b<std::vector<double> > stdSolver(5);
    stdSolver.setVerboseLevel(-1);

    l_bfgs_b<Eigen::VectorXd> eigenSolver(5);
    eigenSolver.setVerboseLevel(-1);
    Eigen::VectorXd eigenX(2), lb(2), ub(2);

    std::cout << "*** MINIMIZING GOLDSTEIN (with armadillo) ***" << std::endl;
    goldstein_price_function< arma::vec > gold;
    gold.set_lower_bound({-2, -2});
    gold.set_upper_bound({2, 2});
    armaX = {-1, 1.5};
    armaSolver.optimize(gold, armaX);
    std::cout << "MINIMUM LOCATED AT : (" << armaX[0] << ", " << armaX[1] << ")" << std::endl << std::endl;

    std::cout << "*** MINIMIZING BOOTH (with std::vector) ***" << std::endl;
    booth_function<std::vector<double> > bl;
    bl.set_lower_bound({-2, -2});
    bl.set_upper_bound({5, 4});
    stdX = {-1, 2};
    stdSolver.optimize(bl, stdX);
    std::cout << "MINIMUM LOCATED AT : (" << stdX[0] << ", " << stdX[1] << ")" << std::endl << std::endl;

    std::cout << "*** MINIMIZING MATYAS (with Eigen3) ***" << std::endl;
    matyas_function< Eigen::VectorXd > matyas;
    lb << -10, -10;
    ub << 10, 10;
    matyas.set_lower_bound(lb);
    matyas.set_upper_bound(ub);
    eigenX << 8, -8;
    eigenSolver.optimize(matyas, eigenX);
    std::cout << "MINIMUM LOCATED AT : (" << eigenX[0] << ", " << eigenX[1] << ")" << std::endl << std::endl;

    std::cout << "*** MINIMIZING BEALE ***" << std::endl;
    beale_function<std::vector<double> > beale;
    beale.set_lower_bound({-4.5, -4.5});
    beale.set_upper_bound({4.5, 4.5});
    stdX = {-1, 0};
    stdSolver.optimize(beale, stdX);
    std::cout << "MINIMUM LOCATED AT : (" << stdX[0] << ", " << stdX[1] << ")" << std::endl << std::endl;

    std::cout << "*** MINIMIZING ROSENBROCK ***" << std::endl;
    armaSolver.setVerboseLevel(0);
    armaSolver.setMachinePrecisionFactor(1e2);
    rosenbrock_function<arma::vec> rf(3);
    rf.set_lower_bound({-10, -10, -10});
    rf.set_upper_bound({10, 10, 10});
    armaX = {-1, 2, 2};
    armaSolver.optimize(rf, armaX);
    std::cout << "***************************" << std::endl;
    std::cout << "MINIMUM LOCATED AT : (" << armaX[0] << ", " << armaX[1] << ", " << armaX[2] << ")" << std::endl;
    return 0;
}

