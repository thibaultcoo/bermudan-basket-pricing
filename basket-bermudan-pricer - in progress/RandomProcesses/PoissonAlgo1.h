#pragma once
#include "Poisson.h"

// First algorithm presented in class to generate poisson distributions

class PoissonAlgo1 : public Poisson
{

public:
	PoissonAlgo1();
	PoissonAlgo1(double _lambda, UniformGenerator* _gen);
	virtual double Generate();

};

