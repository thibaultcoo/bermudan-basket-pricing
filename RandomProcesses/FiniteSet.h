#pragma once
#include "DiscreteGenerator.h"
#include <vector>

//Distribution over a given finite set

class FiniteSet : public DiscreteGenerator
{

public:
	FiniteSet();
	FiniteSet(UniformGenerator* _gen, std::vector<double> _values, std::vector<double> _probas);
	virtual double Generate();

private:
	std::vector<double> values;
	std::vector<double> probas;

};