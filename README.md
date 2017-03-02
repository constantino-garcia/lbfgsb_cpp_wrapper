# L-BFGS-B wrapper for C++
`lbfgsb_cpp_wrapper` is a simple C++ wrapper around the original [Fortran
L-BGSG-B routine](http://users.iems.northwestern.edu/~nocedal/lbfgsb.html), one
of the most widely-used limited-memory quasi-Newton algorithms for
bound-constrained optimization. It tries to be compatible with all king of modern vector-like containers, like [std::vector](http://en.cppreference.com/w/cpp/container/vector), [std::array](http://en.cppreference.com/w/cpp/container/array), [armadillo](http://arma.sourceforge.net/docs.html) vectors or 
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) vectors.

## A quick example
A simple example may be found under the `examples` folder. We give a short overview of the keys steps required to get your code running.

The files required to use the `l_bfgs_b` wrapper class are located under the `include` directory (C++ files) and the `Lbfgsb.3.0` directory (Fortran files). To specify the problem to be optimized you need to implement a tiny class extending the templated-class `problem`, which can be found under the `include` directory.

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

To run the code, it is necessary to link the Fortran routines with the C++ files. Fortunately, it is possible to do it in a straightforward manner using `cmake`.

```bash
mkdir build
cd build
cmake ..
make 
../bin/simple_example
```



## A full example

A full example using all kind of vector-like containers is provided in `l_bfgs_b_example.cpp`. To run the full example, you will need to download and install the latest versions of [armadillo](http://arma.sourceforge.net/docs.html) and [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), since the example makes use of them. Additionally, you will need to set the environment variable `EIGEN3_INCLUDE_DIR` to wherever you had installed `Eigen`. 

After the installation process, you can compile the example using `cmake`:

```bash
# move to the same build directory as in the previous code chunk ...
cd build
# ... and clean it
rm -rf *
# change the installation path if needed
export EIGEN3_INCLUDE_DIR="/usr/local/include/eigen3/"
cmake -DBUILD_SIMPLE_EX=off -DBUILD_FULL_EX=on ..
make
../bin/full_example 
```

## Install it and compile your code

It is possible to use `cmake` to install the library and the required C++ headers:

```bash
cd build
rm -rf *
cmake ..
# Installation requires root permissions
sudo make install
```

If you don't have root permissions or you want to specify an installation folder you
may use:

```bash
cd build
rm -rf *
cmake -DCMAKE_INSTALL_PREFIX:PATH=/path/to/folder ..
make install
```

After the installation, remember to add the shared-library folder to your `LD_LIBRARY_PATH` environment variable. It is possible now to compile your programs without `cmake`:

```bash
# Compile simple_example.cpp
cd examples
# change the path to your library if needed
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/usr/local/lib/lbfgsb_cpp"
g++ simple_example.cpp -std=c++11 -llbfgsb_cpp -L/usr/local/lib/lbfgsb_cpp/ \
    -o simple_example
./simple_example
```

## License
This project is licensed under the terms of the [MPL2.0](https://www.mozilla.org/en-US/MPL/2.0/) license.