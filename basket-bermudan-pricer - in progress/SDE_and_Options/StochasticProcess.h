#pragma once
#include "RandomGenerator.h"
#include "UniformGenerator.h"
#include "SinglePath.h"

class StochasticProcess
{
public:
	StochasticProcess();
	StochasticProcess(RandomGenerator* _gen, int _dim);
	virtual void Simulate(double start_time, double end_time, size_t nb_steps) = 0;
	virtual void Simulate_Antithetic(double start_time, double end_time, size_t nb_steps) = 0;
	virtual void Simulate_VDC(double start_time, double end_time, size_t nb_steps, myLong sim, myLong nbSim) = 0;
	void add_path(SinglePath* Path);
	SinglePath* GetPath(int dimension = 0);
	const double Get_Value(double time, int dim = 0);
	const double Get_Value_antithetic(double time, int dim = 0);
	const std::vector<double> Get_ValueND(double time);
	const std::vector<double> Get_ValueND_antithetic(double time);

protected:
	std::vector<SinglePath*> paths;
	std::vector<SinglePath*> paths_antithetic;
	RandomGenerator* gen;
	int dim;

};

