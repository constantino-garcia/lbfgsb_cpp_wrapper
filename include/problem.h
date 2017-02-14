#ifndef LBFGSB_ADAPTER_PROBLEM_H
#define LBFGSB_ADAPTER_PROBLEM_H

template <class T>
class problem {
public:
     problem(int inputDimension, const T& lowerBound, const T& upperBound) :
            mInputDimension(inputDimension),
            mLowerBound(lowerBound),
            mUpperBound(upperBound){
    };
    problem(int inputDimension) :
            mInputDimension(inputDimension),
            mLowerBound(inputDimension),
            mUpperBound(inputDimension){
        for (int i = 0; i < mInputDimension; ++i) {
            mLowerBound[i] = -std::numeric_limits<double>::infinity();
            mUpperBound[i] = std::numeric_limits<double>::infinity();
        }
    };
    virtual ~problem() = default;
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

};

#endif //LBFGSB_ADAPTER_PROBLEM_H
