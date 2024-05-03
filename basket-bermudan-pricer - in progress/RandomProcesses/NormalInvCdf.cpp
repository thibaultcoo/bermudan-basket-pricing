#include "pch.h"
#include "NormalInvCdf.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

NormalInvCdf::NormalInvCdf()
{
}

NormalInvCdf::NormalInvCdf(double _mean, double _var, UniformGenerator* _gen)
    :Normal(_mean, _var, _gen)
{
}

double RationalApproximation(double t)
{
    // Abramowitz and Stegun formula
    double c[] = { 2.515517, 0.802853, 0.010328 };
    double d[] = { 1.432788, 0.189269, 0.001308 };
    return t - ((c[2] * t + c[1]) * t + c[0]) /
        (((d[2] * t + d[1]) * t + d[0]) * t + 1.0);
}

double NormalCDFInverse(double p)
{
    if (p <= 0.0 || p >= 1.0)
    {
        std::stringstream os;
        os << "Invalid input argument (" << p
            << "); must be larger than 0 but less than 1.";
        throw std::invalid_argument(os.str());
    }
    if (p < 0.5)
    {
        return -RationalApproximation(sqrt(-2.0 * log(p)));
    }
    else
    {
        return RationalApproximation(sqrt(-2.0 * log(1 - p)));
    }
}

double NormalInvCdf::Generate()
{
    double Uniform = generator->Generate();
    return NormalCDFInverse(Uniform);
}