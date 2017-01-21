#ifndef __RANDOMUNIFORM_HPP__
#define __RANDOMUNIFORM_HPP__

#include <stdlib.h>
#include <math.h>
#include <random>
#include <algorithm>

/*
 * BASIC RANDOM ROUTINES COLLECTED IN ONE FILE FOR SIMPLICITY.
 */


class Random
{
public:
    std::random_device r;
    std::seed_seq seed;
    std::mt19937_64 mersenne;
    std::uniform_real_distribution<double> uniform_dist;
    
    Random()
    {
        std::seed_seq _seed{r(), r(), r(), r(), r(), r(), r(), r()};
        std::mt19937_64 _mersenne(_seed);
        std::uniform_real_distribution<double> _uniform_dist(0., 1.);
        seed = _seed;
        mersenne = _mersenne;
        uniform_dist = _uniform_dist;
    }
    
    double get_random_uniform()
    {
        return uniform_dist(mersenne);
    }
    
    template <typename T> inline void shuffle_array(T* array, unsigned size)
    {
        std::shuffle(array, array + size, mersenne);
    }
    
    inline int random_integer(const int min, const int max)
    {
        std::uniform_int_distribution<> rnd_int_gen(min, max);
        return rnd_int_gen(mersenne);
    }
    
    double get_gaussian()
    {
        
        
    }
    
};

#endif