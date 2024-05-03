#include "pch.h"
#include "ClosedFormulas.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>


double norm_cdf(double x)
{
    return 0.5 * std::erfc(-x / std::sqrt(2));
}


double bs_price_call(double S, double strike, double volatility, double maturity, double r)
{
    double stddev = volatility * std::sqrt(maturity);
    double d1 = (std::log(S / strike) + (r + 0.5 * pow(volatility, 2)) * maturity) / stddev;
    double d2 = d1 - stddev;

    double price = S * norm_cdf(d1) - strike * std::exp(-r * maturity) * norm_cdf(d2);

    return price;

}