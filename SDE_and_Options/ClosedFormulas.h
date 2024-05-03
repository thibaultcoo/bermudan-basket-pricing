#pragma once
#include <eigen-3.4.0/Eigen/Dense>
#include <vector>
// Closed formula for the call price in Black-Scholes framework

double bs_price_call(double S, double strike, double volatility, double maturity, double r);

