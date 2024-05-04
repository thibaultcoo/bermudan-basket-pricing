#pragma once
#include "Options.h"
#include "instantiating.h"
#include <eigen-3.4.0/Eigen/Dense>


class EuropeanBasket : public Options
{
	private:
		std::vector<double> weights;
		std::vector<double> spots;
		Eigen::MatrixXd corrMatrix;

	public:
		// Constructors
		EuropeanBasket(StochasticProcess* process, const Instance& Inst);

		// Pricing methods
		double price(int nbSim);
		double priceAntithetic(int nbSim);
		double priceControlVariate(int nbSim);
		double priceQuasiMC(int nbSim);
};