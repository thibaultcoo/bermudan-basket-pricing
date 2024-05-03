#pragma once
#include "Normal.h"

// Inverse Cdf approximation for Normal distribution simulation, particularly useful when using Van Der Corput sequence

class NormalInvCdf :public Normal
{
public:
	NormalInvCdf();
	NormalInvCdf(double _mean, double _var, UniformGenerator* _gen);
	virtual double Generate();

};

double RationalApproximation(double t);
double NormalCDFInverse(double p);

