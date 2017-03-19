/*
 * Copyright Constantino Antonio Garcia 2017
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

//TODO: unify header guards
#ifndef LBFGSB_CPP_RANDOM_VECTOR_GENERATOR_H
#define LBFGSB_CPP_RANDOM_VECTOR_GENERATOR_H

#include <random>

// Use the Curiosly repeating pattern to avoid code duplication
template<class T, typename derived>
class random_vector_generator_base {
public:
    random_vector_generator_base() = default;

    random_vector_generator_base(int size, double minValue, double maxValue) :
            mSize(size),
            mDis(minValue, maxValue) {
        //Will be used to obtain a seed for the random number engine
        std::random_device rd;
        mGen.seed(rd());
    }

    random_vector_generator_base(int size, double minValue, double maxValue, double seed) :
            mSize(size), mGen(seed),
            mDis(minValue, maxValue) {
    }

    ~random_vector_generator_base() = default;

    T operator()() {
        T x(mSize);
        for (int i = 0; i < mSize; i++) {
            x[i] = mDis(mGen);
        }
        return x;
    }

protected:
    int mSize;
    std::mt19937 mGen;
    std::uniform_real_distribution<> mDis;
};


// general version
template<class T>
class random_vector_generator : public random_vector_generator_base<T, random_vector_generator<T> > {
private:
    typedef random_vector_generator_base<T, random_vector_generator<T> > base;

public:
    random_vector_generator() : base() {
    }

    random_vector_generator(int size, double minValue, double maxValue)
            : base(size, minValue, maxValue) {
    }

    random_vector_generator(int size, double minValue, double maxValue, double seed)
            : base(size, minValue, maxValue, seed) {
    }
};


// specialization for std::array
template<typename U, std::size_t N>
class random_vector_generator<std::array<U, N> > :
        public random_vector_generator_base<std::array<U, N>, random_vector_generator<std::array<U, N> > > {
private:
    typedef random_vector_generator_base<std::array<U, N>, random_vector_generator<std::array<U, N> > > base;

public:
    random_vector_generator() : base() {
    }

    random_vector_generator(int size, double minValue, double maxValue)
            : base(size, minValue, maxValue) {
    }

    random_vector_generator(int size, double minValue, double maxValue, double seed)
            : base(size, minValue, maxValue, seed) {
    }

    std::array<U, N> operator()() {
        assert(N == this->mSize);
        std::array<U, N> x;
        std::generate(x.begin(), x.end(),  [&](){return this->mDis(this->mGen);});
        return x;
    }
};


#endif //LBFGSB_CPP_RANDOM_VECTOR_GENERATOR_H
