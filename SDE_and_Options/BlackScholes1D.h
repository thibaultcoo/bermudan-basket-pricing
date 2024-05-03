#pragma once
#include "StochasticProcess.h"

//Simulation in dimension 1

class BlackScholes1D : public StochasticProcess
{
public:
	BlackScholes1D();
	BlackScholes1D(RandomGenerator* _gen, double _s, double _r, double _vol);

protected:
	double s;
	double r;
	double vol;
};