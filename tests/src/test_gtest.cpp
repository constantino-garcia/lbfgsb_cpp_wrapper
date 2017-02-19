#include "gtest/gtest.h"
#include "test_functions.h"
#include "../../include/l_bfgs_b.h"
#include <vector>
#include <armadillo>
#include <Eigen/Dense>
#include <initializer_list>
#include <array>

// Define a test fixture class template that creates the proper
// solver and auxiliar containers for the initial point, the solution and
// the bounds
template<class T>
class l_bfgs_b_test : public testing::Test {
protected:
    // set solvers to work with high-accuracy
    l_bfgs_b_test() : mSolver(l_bfgs_b<T>(5, 5000, 1e1, 1e-12) ){}

    ~l_bfgs_b_test() = default;

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

    l_bfgs_b<T> mSolver;
};


template<class U, std::size_t N>
class l_bfgs_b_test< std::array<U,N> > : public testing::Test {
protected:
    l_bfgs_b_test() : mSolver(l_bfgs_b< std::array<U,N> >()) {}

    ~l_bfgs_b_test() = default;

    std::array<U, N> fill_proper_container(const std::initializer_list<double>& initList) {
        assert(N == initList.size());
        std::array<U, N> x;
        std::copy(initList.begin(), initList.end(), std::begin(x));
        return x;
    }

    l_bfgs_b< std::array<U,N> > mSolver;
};


using testing::Types;
typedef Types<std::vector<double>, arma::vec, Eigen::VectorXd, std::array<double, 2> > Implementations;


TYPED_TEST_CASE(l_bfgs_b_test, Implementations);

TYPED_TEST(l_bfgs_b_test, rosenbrock) {
    int n = 2;
    rosenbrock_function<TypeParam> func(n);
    TypeParam x = this->fill_proper_container({-1,2});
    TypeParam sol = this->fill_proper_container({1,1});
    TypeParam lb = this->fill_proper_container({-10,-10});
    TypeParam ub = this->fill_proper_container({10,10});
    func.set_lower_bound(lb);
    func.set_upper_bound(ub);
    this->mSolver.optimize(func, x);
    for (int i = 0; i < sol.size(); ++i) {
        ASSERT_NEAR(sol[i], x[i], 1e-5);
    }
}

TYPED_TEST(l_bfgs_b_test, beale) {
    beale_function<TypeParam> func;
    TypeParam x = this->fill_proper_container({-1,0});
    TypeParam sol = this->fill_proper_container({3,0.5});
    TypeParam lb = this->fill_proper_container({-4.5,-4.5});
    TypeParam ub = this->fill_proper_container({4.5,4.5});
    func.set_lower_bound(lb);
    func.set_upper_bound(ub);
    this->mSolver.optimize(func, x);
    for (int i = 0; i < sol.size(); ++i) {
        ASSERT_NEAR(sol[i], x[i], 1e-5);
    }
}


TYPED_TEST(l_bfgs_b_test, booth) {
    booth_function<TypeParam> func;
    TypeParam x = this->fill_proper_container({8,9.7});
    TypeParam sol = this->fill_proper_container({1,3});
    TypeParam lb = this->fill_proper_container({-10,-10});
    TypeParam ub = this->fill_proper_container({10,10});
    func.set_lower_bound(lb);
    func.set_upper_bound(ub);
    this->mSolver.optimize(func, x);
    for (int i = 0; i < sol.size(); ++i) {
        ASSERT_NEAR(sol[i], x[i], 1e-5);
    }
}



TYPED_TEST(l_bfgs_b_test, matyas) {
    matyas_function<TypeParam> func;
    TypeParam x = this->fill_proper_container({8,-9.7});
    TypeParam sol = this->fill_proper_container({0,0});
    TypeParam lb = this->fill_proper_container({-10,-10});
    TypeParam ub = this->fill_proper_container({10,10});
    func.set_lower_bound(lb);
    func.set_upper_bound(ub);
    this->mSolver.optimize(func, x);
    for (int i = 0; i < sol.size(); ++i) {
        ASSERT_NEAR(sol[i], x[i], 1e-5);
    }
}
