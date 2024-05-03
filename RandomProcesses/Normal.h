#pragma once
#include "ContinuousGenerator.h"

// Normal distribution generations

class Normal : public ContinuousGenerator
{

public:
	Normal();
	Normal(double _mean, double _var, UniformGenerator* _gen);

protected:
	double mean;
	double var;

};

