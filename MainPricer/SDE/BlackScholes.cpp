#include "pch.h"
#include "BlackScholes.h"

BlackScholes::BlackScholes() { }

BlackScholes::BlackScholes(RandomGenerator* generator, std::vector<double> spots, std::vector<double> rates, Eigen::MatrixXd corrMatrix, int dimension)
	: StochasticProcess(generator, dimension), spots(spots), rates(rates), corrMatrix(corrMatrix) { }