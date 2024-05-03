#pragma once
#include "Poisson.h"
#include "Exponential.h"

//Second method presented in class to simulate Poisson distribution using exponential distribution

class PoissonAlgo2 : public Poisson
{

public:
	PoissonAlgo2();
	PoissonAlgo2(Exponential* _gen);
	virtual double Generate();

private:
	Exponential* generator;


};

