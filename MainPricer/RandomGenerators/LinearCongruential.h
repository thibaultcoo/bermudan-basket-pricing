#pragma once
#include "PseudoGenerator.h"

class LinearCongruential : public PseudoGenerator
{
	protected:
		myLong Multiplier;
		myLong Increment;
		myLong Modulus;

	public:
		LinearCongruential();
		LinearCongruential(myLong _seed, myLong _multiplier, myLong _increment, myLong _modulus);
		double Generate();
		myLong get_Modulus();
};