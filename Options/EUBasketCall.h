#pragma once
#include "Options.h"
#include <eigen-3.4.0/Eigen/Dense>


class EUBasketCall : public Options
{
	private:
		std::vector<double> weights;
		std::vector<double> spots;
		Eigen::MatrixXd corrMatrix;

	public:
		EUBasketCall();
		EUBasketCall(StochasticProcess* process, double strike, std::vector<double> rates, double maturity, std::vector<double> weights);
		EUBasketCall(StochasticProcess* process, double strike, std::vector<double> rates, double maturity, std::vector<double> weights, std::vector<double> spots, Eigen::MatrixXd corrMatrix);
		double price(int nbSim);
		double priceAntithetic(int nbSim);
		double priceControlVariate(int nbSim);
		double priceVDC(int nbSim);
};

double norm_cdf(double x);
double BS_Call(double spot, double strike, double volatility, double maturity, double rate);
double compute_expected_value_control_variate(std::vector<double> spots, std::vector<double> weights, double strike, double rate, Eigen::MatrixXd corrMatrix, double maturity);
