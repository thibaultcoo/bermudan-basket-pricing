#pragma once
#include "Options.h"
#include "instantiating.h"
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
		// Constructor
		BermudanBasket(StochasticProcess* process, const Instance& Inst, int L = 5);

		// Pricing methods
		double price(int NbSim);
		double priceAntithetic(int NbSim);
		double priceControlVariate(int NbSim);
		double priceQuasiMC(int NbSim);

		// Helper Functions
		void backwardInduction(int dateIndex, int numSimulations, Eigen::MatrixXd& basketValues, Eigen::MatrixXd& exerciseTimes);
		double discountedPayoff(int simIndex, int numExerciseDates, Eigen::MatrixXd& basketValues, Eigen::MatrixXd& exerciseTimes);
};