#pragma once
#include "PseudoGenerator.h"
#include "LinearCongruential.h"


class EcuyerCombined : public PseudoGenerator
{
	protected:
		LinearCongruential Generator1;
		LinearCongruential Generator2;

	public:
		EcuyerCombined();
		double Generate();
};