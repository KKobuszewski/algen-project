#ifndef __MATH_EXTENSIONS_H__
#define __MATH_EXTENSIONS_H__


#include <math.h>

inline unsigned int factorial(const unsigned int n)
{
    unsigned int retval = 1;
    for (int i = n; i > 1; --i) retval *= i;
    return retval;
}


#endif