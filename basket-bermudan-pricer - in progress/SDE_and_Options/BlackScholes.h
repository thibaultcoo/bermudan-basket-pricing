#pragma once
#include "StochasticProcess.h"
#include <vector>
#include <iostream>
#include <eigen-3.4.0/Eigen/Dense>

//Simulations in dimension n

class BlackScholes : public StochasticProcess
{
public:
	BlackScholes();
	BlackScholes(RandomGenerator* _gen, std::vector<double> _s, std::vector<double> _r, Eigen::MatrixXd _VarCovar, int dim);

protected:
	std::vector<double> s;
	std::vector<double> r;
	Eigen::MatrixXd VarCovar; //Variance covariance matrix
};