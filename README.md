# L-BFGS-B wrapper for C++
`lbfgsb_cpp_wrapper` is a simple C++ wrapper around the original [Fortran
L-BGSG-B routine](http://users.iems.northwestern.edu/~nocedal/lbfgsb.html), one
of the most widely-used limited-memory quasi-Newton algorithms for
bound-constrained optimization. It tries to be compatible with all king of modern vector-like containers, like [std::vector](http://en.cppreference.com/w/cpp/container/vector), [std::array](http://en.cppreference.com/w/cpp/container/array), [armadillo](http://arma.sourceforge.net/docs.html) vectors or 
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) vectors.

## A quick example
A simple example may be found under the `simple_example` folder. We give a short overview of the keys steps required to get your code running.

The files required to use the `l_bfgs_b` wrapper class are located under the *include* directory (C++ files) and the `Lbfgsb.3.0` directory (Fortran files). To specify the problem to be optimized you need to implement a tiny class extending the templated-class `problem`, which can be found under the `include` directory.

```c++
template<class T>
class quadratic_problem : public problem<T> {
public:
    // here, 2 is the problem input dimension
    quadratic_problem() : problem(2){}

    double operator()(const T& x) {
        return std::pow(x[0] - 0.5, 2) +  std::pow(x[1] - 1, 2);
    }

    void gradient(const T& x, T& gr) {
        gr[0] = 2 * (x[0] - 0.5);
        gr[1] = 2 * (x[1] - 1);
    }
};
```

Note that this class specifies the dimension of the problem through the constructor, 
the objective function and its gradient. To specify the box-constraints of the problems we can use the `set_lower_bound` and `set_upper_bound` methods inherited from class `problem`:

```c++
 typedef std::array<double,2> my_vector;
 quadratic_problem<my_vector> qp;
 qp.set_lower_bound({-2, -2});
 qp.set_upper_bound({2, 2});
```

The optimization is then performed using the `l_bfgs_b` class:

```c++
 my_vector initPoint = {2, 3};
 l_bfgs_b<my_vector> solver;
 solver.optimize(qp, initPoint);
 std::cout << "MINIMUM LOCATED AT : (" << 
           initPoint[0] << ", " << initPoint[1] << ")" << 
           std::endl;

```

To run the code, it is necessary to link the Fortran files with the C++ files. Fortunately, it is possible to do it in a straightforward manner using `cmake`.

```bash
cd simple_example
cmake CMakeLists.txt 
make 
./simple_example
```



## A full example

A full example using all kind of vector-like containers is provided in `l_bfgs_b_example.cpp`. To run the full example, you will need to download and install the latest versions of [armadillo](http://arma.sourceforge.net/docs.html) and [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), since the example makes use of them. Additionally, you will need to set the environment variable `EIGEN3_INCLUDE_DIR` to wherever you had installed `Eigen`. 

After the installation process, you can compile the example using `cmake`:

```bash
# change the installation path if needed
export EIGEN3_INCLUDE_DIR="/usr/local/include/eigen3/"
cmake CMakeLists.txt 
make example
```
