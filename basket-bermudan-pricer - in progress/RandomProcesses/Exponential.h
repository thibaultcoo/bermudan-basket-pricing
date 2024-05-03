#pragma once
#include "ContinuousGenerator.h" 

//Exponential distribution

class Exponential :public ContinuousGenerator
{
public:
	Exponential();
	Exponential(double _lambda, UniformGenerator* _gen);

protected:
	double lambda;
};

