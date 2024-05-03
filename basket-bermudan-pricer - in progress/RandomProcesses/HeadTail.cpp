#include "pch.h"
#include "HeadTail.h"
#include "LinearCongruential.h"

HeadTail::HeadTail()
{
}

HeadTail::HeadTail(UniformGenerator* _gen)
	: DiscreteGenerator(_gen)
{
}

double HeadTail::Generate()
{
	double u = generator->Generate();
	if (u <= 0.5)
		return 1.;
	else
		return 0.;
}