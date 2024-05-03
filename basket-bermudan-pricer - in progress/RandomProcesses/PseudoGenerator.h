#pragma once
#include "UniformGenerator.h"

class PseudoGenerator : public UniformGenerator
{
public:
	PseudoGenerator();
	PseudoGenerator(myLong inputSeed);

protected:
	myLong seed;
	myLong currentNumber;
};

