#include "pch.h"
#include "Binomial.h"

Binomial::Binomial(UniformGenerator* _gen, myLong _n, double _p)
	: DiscreteGenerator(_gen), n(_n), p(_p)
{

}

double Binomial::Generate()
{
	double u;
	double result = 0;

	for (myLong i = 0; i < n; ++i)
	{
		u = generator->Generate();
		if (u <= p)
			result += 1.;
	}
	return result;
}