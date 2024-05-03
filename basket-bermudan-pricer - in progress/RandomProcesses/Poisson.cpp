#include "pch.h"
#include "Poisson.h"

Poisson::Poisson()
{
}

Poisson::Poisson(double _lambda, UniformGenerator* _gen)
	:DiscreteGenerator(_gen), lambda(_lambda)
{
}