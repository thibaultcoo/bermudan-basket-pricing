#pragma once
#include "Normal.h"

// Use of Central Limit Theorem to get Normal distribution

class NormalCLT : public Normal
{
public:
	NormalCLT();
	NormalCLT(double _mean, double _var, UniformGenerator* _gen);
	virtual double Generate();
};

