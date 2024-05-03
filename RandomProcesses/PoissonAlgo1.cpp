#include "pch.h"
#include "PoissonAlgo1.h"
#include <math.h>
#include "iostream"

PoissonAlgo1::PoissonAlgo1()
{
}

PoissonAlgo1::PoissonAlgo1(double _lambda, UniformGenerator* _gen)
	:Poisson(_lambda, _gen)
{
}

double PoissonAlgo1::Generate()
{
	double u = generator->Generate();

	int k = 0;
	myLong factorial = 1;
	double p;
	double P = exp(-1 * lambda) * pow(lambda, k) / factorial;

	while ((P < u))
	{
		p = exp(-1 * lambda) * pow(lambda, k) / factorial;
		P += lambda * p / (k + 1);
		k++;
		factorial *= k;
	}

	return k;

}