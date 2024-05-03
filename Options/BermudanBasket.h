#pragma once
#include "Options.h"
#include <eigen-3.4.0/Eigen/Dense>


class BermudanBasket : public Options
{
	private:
		std::vector<double> exerciseDates;
		int L;
		std::vector<double> weights;
		std::vector<double> spots;
		Eigen::MatrixXd corrMatrix;

	public:
		BermudanBasket();
		BermudanBasket(StochasticProcess* process, double strike, std::vector<double> rates, double maturity, std::vector<double> weights, std::vector<double> exerciseDates, int L = 5);
		BermudanBasket(StochasticProcess* process, double strike, std::vector<double> rates, double maturity, std::vector<double> weights, std::vector<double> exerciseDates, std::vector<double> spots, Eigen::MatrixXd corrMatrix, int L = 5);
		double price(int NbSim);
		double priceAntithetic(int NbSim);
		double priceControlVariate(int NbSim);
		double priceVDC(int NbSim);
};

double Compute_E_Ybis(std::vector<double> S, std::vector<double> weights, double K, double r, Eigen::MatrixXd VarCovar, double T);
double norm_cdf(double x);
double BS_Call(double spot, double strike, double volatility, double maturity, double rate);

