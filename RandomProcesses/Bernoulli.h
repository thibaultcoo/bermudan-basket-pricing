#pragma once
#include "DiscreteGenerator.h"
//Bernoulli law with parameter p, p€[0,1]

class Bernoulli : public DiscreteGenerator
{
public:
	Bernoulli();
	Bernoulli(UniformGenerator* _gen, double _p);
	virtual double Generate();

private:
	double p;
};

