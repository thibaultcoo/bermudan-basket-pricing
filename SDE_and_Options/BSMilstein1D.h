#pragma once
#include "BlackScholes1D.h"

// Simulation using the Milstein scheme in one dimension

class BSMilstein1D : public BlackScholes1D
{
public:
	BSMilstein1D();
	BSMilstein1D(RandomGenerator* _gen, double _s, double _r, double _vol);
	void Simulate(double start_time, double end_time, size_t nb_steps);
	void Simulate_Antithetic(double start_time, double end_time, size_t nb_steps);
	void Simulate_VDC(double start_time, double end_time, size_t nb_steps, myLong sim, myLong nbSim);

private:
	std::vector<double> v_VDC;
};