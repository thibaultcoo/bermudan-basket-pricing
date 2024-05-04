#pragma once
#include "BlackScholes.h"

class BSMilstein : public BlackScholes
{
	public:
		BSMilstein();
		BSMilstein(RandomGenerator* generator, std::vector<double> spots, std::vector<double> rates, Eigen::MatrixXd corrMatrix, int dimension);
		void Simulate(double start_time, double end_time, size_t nb_steps);
		void SimulateAntithetic(double start_time, double end_time, size_t nb_steps);
		void SimulateQuasiMC(double start_time, double end_time, size_t nb_steps, myLong sim, myLong nbSim);

	private:
		Eigen::MatrixXd choleskyDecomposition;
		Eigen::MatrixXd vanDerCorputMatrix;
};