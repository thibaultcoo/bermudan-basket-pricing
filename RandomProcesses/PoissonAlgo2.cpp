#include "pch.h"
#include "PoissonAlgo2.h"


PoissonAlgo2::PoissonAlgo2()
{
}

PoissonAlgo2::PoissonAlgo2(Exponential* _gen)
{
	generator = _gen;
}

double PoissonAlgo2::Generate()
{
	double sum = 0;
	myLong k = 0;

	while (sum < 1)
	{
		sum += generator->Generate();
		k++;
	}

	return k - 1;

}