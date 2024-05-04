#pragma once
#include "StochasticProcess.h"
#include <eigen-3.4.0/Eigen/Dense>
#include <vector>

class BlackScholes : public StochasticProcess
{
	public:
		BlackScholes();
		BlackScholes(RandomGenerator* generator, std::vector<double> spots, std::vector<double> rates, Eigen::MatrixXd corrMatrix, int dimension);

	protected:
		std::vector<double> spots;
		std::vector<double> rates;
		Eigen::MatrixXd corrMatrix;
};