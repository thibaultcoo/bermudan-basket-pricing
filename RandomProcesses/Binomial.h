#pragma once
#include "DiscreteGenerator.h"

// Binomial law of parameter n € N (integer) and p€[0,1]

class Binomial : public DiscreteGenerator
{
public:
	Binomial();
	Binomial(UniformGenerator* _gen, myLong _n, double _p);
	virtual double Generate();

private:
	double p;
	myLong n;
};