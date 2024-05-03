#pragma once
#include "RandomGenerator.h"
#include "UniformGenerator.h"

class ContinuousGenerator : public RandomGenerator
{
public:
	ContinuousGenerator();
	ContinuousGenerator(UniformGenerator* _gen);

protected:
	UniformGenerator* generator;
};

