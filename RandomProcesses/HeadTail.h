#pragma once
#include "DiscreteGenerator.h"

//Head tail distribution

class HeadTail : public DiscreteGenerator
{
public:
	HeadTail();
	HeadTail(UniformGenerator* _gen);
	virtual double Generate();

};