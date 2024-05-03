#include "pch.h"
#include "Bernoulli.h"

Bernoulli::Bernoulli(UniformGenerator* _gen, double _p)
	:DiscreteGenerator(_gen), p(_p)
{

}

double Bernoulli::Generate()
{
	double u = generator->Generate();
	if (u <= p)
		return 1.;
	else
		return 0.;
}