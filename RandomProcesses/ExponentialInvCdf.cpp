#include "pch.h"
#include "ExponentialInvCdf.h"
#include <math.h>

ExponentialInvCdf::ExponentialInvCdf()
{
}

ExponentialInvCdf::ExponentialInvCdf(double _lambda, UniformGenerator* _gen)
	:Exponential(_lambda, _gen)
{
}

double ExponentialInvCdf::Generate()
{
	double u = generator->Generate();

	return -log(u) / lambda;
}