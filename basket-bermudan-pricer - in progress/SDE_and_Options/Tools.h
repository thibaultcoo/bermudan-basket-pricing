#pragma once
#include <eigen-3.4.0/Eigen/Dense>
#include <vector>

double norm_cdf(double x);
double fp(double x, double p);
double BS_Call(double spot, double strike, double volatility, double maturity, double rate);
double compute_expected_value_control_variate(std::vector<double> spots, std::vector<double> weights, double strike, double rate, Eigen::MatrixXd corrMatrix, double maturity);