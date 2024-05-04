#include "pch.h"
#include "PseudoGenerator.h"

PseudoGenerator::PseudoGenerator() : seed(0.)
{
	currentNumber = seed;
}

PseudoGenerator::PseudoGenerator(myLong inputSeed) : seed(inputSeed)
{
	currentNumber = inputSeed;
}