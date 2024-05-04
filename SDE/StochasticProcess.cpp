#include "pch.h"
#include "StochasticProcess.h"
#include <vector>
#include "SingleTrajectory.h"

StochasticProcess::StochasticProcess() {}

// Diffuses a stochastic process for the multi-dimensional underlying following Black-Scholes theory
StochasticProcess::StochasticProcess(RandomGenerator* _gen, int _dim)
	: gen(_gen), dim(_dim) {
	trajectory.resize(_dim);
	trajectoryAntithetic.resize(_dim);
	for (int i = 0; i < _dim; ++i) {
		trajectory[i] = new SingleTrajectory();
		trajectoryAntithetic[i] = new SingleTrajectory();
	}
}

// Helper function to append a trajectory to another
void StochasticProcess::addTrajectory(SingleTrajectory* Path) {
	trajectory.push_back(Path);
}

// Helper function to reach to a specific trajectory among a group of trajectories
SingleTrajectory* StochasticProcess::getTrajectory(int dimension) {
	return trajectory[dimension];
}

// Helper function to isolate a specific needed value
const double StochasticProcess::getValue(double time, int dim) {
	return trajectory[dim]->GetValue(time);
}

const double StochasticProcess::getValueAntithetic(double time, int dim) {
	return trajectoryAntithetic[dim]->GetValue(time);
}


const std::vector<double> StochasticProcess::getValueMulti(double time) {
	std::vector<double> res;
	res.reserve(dim);

	for (int i = 0; i < dim; i++) {
		res.push_back(getValue(time, i));
	}
	return res;
}

const std::vector<double> StochasticProcess::getValueAntitheticMulti(double time) {
	std::vector<double> res;
	res.reserve(dim);

	for (int i = 0; i < dim; i++) {
		res.push_back(getValueAntithetic(time, i));
	}
	return res;
}