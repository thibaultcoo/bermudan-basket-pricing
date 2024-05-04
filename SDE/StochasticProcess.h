#pragma once
#include "RandomGenerator.h"
#include "UniformGenerator.h"
#include "SingleTrajectory.h"

class StochasticProcess
{
	public:
		StochasticProcess();
		StochasticProcess(RandomGenerator* _gen, int _dim);

		virtual void Simulate(double start, double end, size_t steps) = 0;
		virtual void SimulateAntithetic(double start, double end, size_t steps) = 0;
		virtual void SimulateQuasiMC(double start, double end, size_t steps, myLong sim, myLong nbSim) = 0;

		void addTrajectory(SingleTrajectory* Path);
		SingleTrajectory* getTrajectory(int dimension = 0);
		const double getValue(double time, int dim = 0);
		const double getValueAntithetic(double time, int dim = 0);

		const std::vector<double> getValueMulti(double time);
		const std::vector<double> getValueAntitheticMulti(double time);

	protected:
		std::vector<SingleTrajectory*> trajectory;
		std::vector<SingleTrajectory*> trajectoryAntithetic;

		RandomGenerator* gen;
		int dim;
};