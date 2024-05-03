#pragma once
#include "DiscreteGenerator.h"

// Poisson distribution with parameter lambda>0

class Poisson : public DiscreteGenerator
{
public:
	Poisson();
	Poisson(double _lambda, UniformGenerator* _gen);

protected:
	double lambda;
};

