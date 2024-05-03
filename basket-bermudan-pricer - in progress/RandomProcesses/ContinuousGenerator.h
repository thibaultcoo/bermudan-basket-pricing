#pragma once
#include "RandomGenerator.h"
#include "UniformGenerator.h"

class ContinuousGenerator : public RandomGenerator
{
	protected:
		UniformGenerator* generator;

	public:
		ContinuousGenerator();
		ContinuousGenerator(UniformGenerator* _gen);
};

