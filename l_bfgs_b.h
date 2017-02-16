#ifndef L_BFGS_B_CPP_WRAPPER_H
#define L_BFGS_B_CPP_WRAPPER_H

#include <cassert>
#include "include/problem.h"

extern "C" {
void setulb_wrapper(int *n, int *m, double x[], double l[], double u[], int nbd[], double *f,
                    double g[], double *factr, double *pgtol, double wa[], int iwa[], int *itask,
                    int *iprint, int *icsave, bool* lsave0, bool* lsave1, bool* lsave2, bool* lsave3,
                    int isave[], double dsave[]);
}

template<class T>
class l_bfgs_b {
public:
    l_bfgs_b(): l_bfgs_b(5) {

    };

    l_bfgs_b(int memorySize): l_bfgs_b(memorySize, 500, 1e7, 1e-15) {

    }

    // Typical values for machinePrecisionFactor : 1e+12 for
    // low accuracy; 1e+7 for moderate accuracy; 1e+1 for extremely
    // high accuracy.
    l_bfgs_b(int memorySize, int maximumNumberOfIterations,
             double machinePrecisionFactor, double projectedGradientTolerance)
    : mMemorySize(memorySize),
            mMaximumNumberOfIterations(maximumNumberOfIterations),
            mMachinePrecisionFactor(machinePrecisionFactor),
            mProjectedGradientTolerance(projectedGradientTolerance),
            mVerboseLevel(-1) {
    }

    ~l_bfgs_b() = default;

    int getMemorySize() const {
        return mMemorySize;
    }

    void setMemorySize(int memorySize) {
        // TODO: check if 1 is a valid input
        if (memorySize <= 1) {
            throw std::invalid_argument("memorySize should be > 1");
        }
        mMemorySize = memorySize;
    }

    double getMachinePrecisionFactor() const {
        return mMachinePrecisionFactor;
    }

    void setMachinePrecisionFactor(double machinePrecisionFactor) {
        if (machinePrecisionFactor <= 0) {
            throw std::invalid_argument("machinePrecisionFactor should be > 0");
        }
        mMachinePrecisionFactor = machinePrecisionFactor;
    }

    double getProjectedGradientTolerance() const {
        return mProjectedGradientTolerance;
    }

    void setProjectedGradientTolerance(double projectedGradientTolerance) {
        mProjectedGradientTolerance = projectedGradientTolerance;
    }

    int getVerboseLevel() const {
        return mVerboseLevel;
    }

    void setVerboseLevel(int verboseLevel) {
        // TODO: discretize verboseLevel as enum class
        mVerboseLevel = verboseLevel;
    }

    void optimize(problem<T>& pb, T &x0) {
        int n = pb.getInputDimension();
        // prepare variables for the algorithm
        std::vector<double> mX(n);
        std::vector<double> mGradient(n);
        std::vector<double> mLowerBound(n);
        std::vector<double> mUpperBound(n);
        std::vector<int> mNbd(n);
        std::vector<double> mWorkArray(2 * mMemorySize * n + 5 * n +
                                       11 * mMemorySize * mMemorySize + 8 * mMemorySize);
        std::vector<int> mIntWorkArray(3 * n);

        T lowerBound = pb.get_lower_bound();
        T upperBound = pb.get_upper_bound();
        bool hasLowerBound, hasUpperBound;
        for (int i = 0; i < n; ++i) {
            mX[i] = x0[i];
            mLowerBound[i] = lowerBound[i];
            mUpperBound[i] = upperBound[i];
            hasLowerBound = !std::isinf(mLowerBound[i]);
            hasUpperBound = !std::isinf(mUpperBound[i]);
            // nbd(i)=0 if x(i) is unbounded,
            // 1 if x(i) has only a lower bound,
            // 2 if x(i) has both lower and upper bounds, and
            // 3 if x(i) has only an upper bound.
            if (hasLowerBound) {
                if (hasUpperBound) {
                    mNbd[i] = 2;
                } else {
                    mNbd[i] = 1;
                }
            } else if (hasUpperBound) {
                mNbd[i] = 3;
            } else {
                mNbd[i] = 0;
            }
        }

        double f = pb(x0);
        T gr;
        pb.gradient(x0, gr);
        fill_pointer(mGradient, gr, n);

        int i = 0;
        int itask = 0;
        int icsave = 0;

        bool test = false;

        while ((i < mMaximumNumberOfIterations) && (
                (itask == 0) || (itask == 1) || (itask == 2) || (itask == 3)
        )) {
            setulb_wrapper(&n, &mMemorySize, &mX[0], &mLowerBound[0], &mUpperBound[0], &mNbd[0], &f,
                           &mGradient[0],
                           &mMachinePrecisionFactor, &mProjectedGradientTolerance,
                           &mWorkArray[0], &mIntWorkArray[0], &itask, &mVerboseLevel,
                           &icsave, &mBoolInformation[0], &mBoolInformation[1],
                           &mBoolInformation[2], &mBoolInformation[3],
                           &mIntInformation[0], &mDoubleInformation[0]);
            // assert that impossible values do not occur
            assert(icsave <= 14 && icsave >= 0);
            assert(itask <= 12 && itask >= 0);

            fill_container(mX, x0, n);
            
            if (itask == 2 || itask == 3) {
                f = pb(x0);
                pb.gradient(x0, gr);
                fill_pointer(mGradient, gr, n);
            }
            
            i = mIntInformation[29];
        }

    }


private:
    int mMemorySize;
    double mMachinePrecisionFactor;
    double mProjectedGradientTolerance;
    int mVerboseLevel;
    int mMaximumNumberOfIterations;
    // interface to Fortran code
    char mTask[60];
    char mCharInformation[60];
    // Use ints to model boolean information to avoid problems with the Fortran-C++ binding
    bool mBoolInformation[4];
    int mIntInformation[44];
    double mDoubleInformation[29];

    // TODO: remove the fill_XXX functions
    // the array and the container should have the proper dimensions
    // there is no resizing inside these functions for efficiency
    void fill_container(const std::vector<double>& in, T &out, int n) {
        for (int i = 0; i < n; ++i) {
            out[i] = in[i];
        }
    }
    void fill_pointer(std::vector<double>& out, const T &in, int n) {
        for (int i = 0; i < n; ++i) {
            out[i] = in[i];
        }
    }
};


#endif //L_BFGS_B_CPP_WRAPPER_H
