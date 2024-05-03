#pragma once
#include "Options.h"
#include <eigen-3.4.0/Eigen/Dense>
#include <vector>


class BermudanBasket : public Options
{
	private:
		std::vector<double> exerciseDates;
		int L;
		std::vector<double> weights;
		std::vector<double> spots;
		Eigen::MatrixXd corrMatrix;

	public:
		// Constructors
		BermudanBasket();
		BermudanBasket(StochasticProcess* process, double strike, std::vector<double> rates, double maturity, std::vector<double> weights, std::vector<double> exerciseDates, int L = 5);
		BermudanBasket(StochasticProcess* process, double strike, std::vector<double> rates, double maturity, std::vector<double> weights, std::vector<double> exerciseDates, std::vector<double> spots, Eigen::MatrixXd corrMatrix, int L = 5);
		
		// Pricing methods
		double price(int NbSim);
		double priceAntithetic(int NbSim);
		double priceControlVariate(int NbSim);
		double priceVDC(int NbSim);

		// Helper Functions
		void backwardInduction(int dateIndex, int numSimulations, Eigen::MatrixXd& basketValues, Eigen::MatrixXd& exerciseTimes);
		double discountedPayoff(int simIndex, int numExerciseDates, Eigen::MatrixXd& basketValues, Eigen::MatrixXd& exerciseTimes);
};