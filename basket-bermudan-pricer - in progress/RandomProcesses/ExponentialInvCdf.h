#pragma once
#include "Exponential.h"

//Using the inverse of the Cdf to simulate an exponential process

class ExponentialInvCdf : public Exponential
{

public:
	ExponentialInvCdf();
	ExponentialInvCdf(double _lambda, UniformGenerator* _gen);
	virtual double Generate();

};