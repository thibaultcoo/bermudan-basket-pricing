#pragma once
#include "UniformGenerator.h"

class PseudoGenerator : public UniformGenerator
{
	protected:
		myLong seed;
		myLong currentNumber;

	public:
		PseudoGenerator();
		PseudoGenerator(myLong inputSeed);
};

