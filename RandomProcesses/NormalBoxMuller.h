#pragma once
#include "Normal.h"

// Use of Box-Muller method to simulate the Normal Distribution

class NormalBoxMuller : public Normal
{

public:
	NormalBoxMuller();
	NormalBoxMuller(double _mean, double _var, UniformGenerator* _gen);
	virtual double Generate();

private:
	bool NewSimulation = true;
	double SecondNormal;
};