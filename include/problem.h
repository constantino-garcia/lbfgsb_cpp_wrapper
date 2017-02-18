#ifndef LBFGSB_ADAPTER_PROBLEM_H
#define LBFGSB_ADAPTER_PROBLEM_H

#include <array>
// TODO: check that inputDimensions matches lowerBound.size()
template <typename T, typename derived>
class problem_base {
public:
     problem_base(int inputDimension, const T& lowerBound, const T& upperBound) :
            mInputDimension(inputDimension),
            mLowerBound(lowerBound),
            mUpperBound(upperBound){
    };
    problem_base(int inputDimension) :
            mInputDimension(inputDimension),
            mLowerBound(inputDimension),
            mUpperBound(inputDimension){
        set_default_bounds();
    }

    virtual ~problem_base() = default;
    int getInputDimension() const {
        return mInputDimension;
    }

    T get_lower_bound() const {
        return mLowerBound;
    }

    void set_lower_bound(const T& lowerBound) {
        mLowerBound = lowerBound;
    }

    T get_upper_bound() const {
        return mUpperBound;
    }

    void set_upper_bound(const T& upperBound) {
        mUpperBound = upperBound;
    }

    virtual double operator()(const T& x) = 0;
    virtual void gradient(const T& x, T& gr) = 0;

protected:
    int mInputDimension;
    T mLowerBound;
    T mUpperBound;

    problem_base() = default;

    void set_default_bounds() {
        for (int i = 0; i < mInputDimension; ++i) {
            mLowerBound[i] = -::std::numeric_limits<double>::infinity();
            mUpperBound[i] = ::std::numeric_limits<double>::infinity();
        }
    };

    static void check_input_dimension(int N, int inputDimension) {
         if (N != inputDimension) {
             throw std::invalid_argument("The container's size does not match the problem's input dimension");
        }
    }
};

// The general version.
template <typename T>
class problem : public problem_base<T, problem<T> > {
private:
    typedef problem_base<T, problem<T> > base;

public:
    problem(int inputDimension, const T& lowerBound, const T& upperBound) :
          base(inputDimension, lowerBound, upperBound) {
    }

    problem(int inputDimension) : base(inputDimension){
    }
};

// Specialization for std::array
template <typename U, std::size_t  N>
class problem< std::array<U, N> > : public problem_base<std::array<U,N>, problem< std::array<U,N> > >{
private:
    typedef problem_base<std::array<U, N>, problem<std::array<U, N> > > base;

public:
    problem(int inputDimension, const std::array<U,N>& lowerBound, const std::array<U,N>& upperBound) :
          base() {
        this->check_input_dimension(N, inputDimension);
        this->mInputDimension = inputDimension;
        this->mLowerBound = lowerBound;
        this->mUpperBound = upperBound;
    }

    problem(int inputDimension) : base(){
        this->check_input_dimension(N, inputDimension);
        this->mInputDimension = inputDimension;
        this->set_default_bounds();
    }
};

#endif //LBFGSB_ADAPTER_PROBLEM_H
