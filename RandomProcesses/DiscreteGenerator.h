#pragma once
#include "UniformGenerator.h"

class DiscreteGenerator : public RandomGenerator
{
public:
	DiscreteGenerator();
	DiscreteGenerator(UniformGenerator* _gen);

protected:
	UniformGenerator* generator;
};

