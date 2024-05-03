#pragma once
#include "Options.h"
#include <eigen-3.4.0/Eigen/Dense>
#include <vector>


class EUBasketCall : public Options
{
	private:
		std::vector<double> weights;
		std::vector<double> spots;
		Eigen::MatrixXd corrMatrix;

	public:
		// Constructors
		EUBasketCall();
		EUBasketCall(StochasticProcess* process, double strike, std::vector<double> rates, double maturity, std::vector<double> weights);
		EUBasketCall(StochasticProcess* process, double strike, std::vector<double> rates, double maturity, std::vector<double> weights, std::vector<double> spots, Eigen::MatrixXd corrMatrix);
		
		// Pricing methods
		double price(int nbSim);
		double priceAntithetic(int nbSim);
		double priceControlVariate(int nbSim);
		double priceVDC(int nbSim);
};