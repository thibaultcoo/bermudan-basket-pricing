#pragma once
#include "QuasiGenerator.h"

//Van Der Corput sequence generation

class VanDerCorput : public QuasiGenerator
{
public:
	VanDerCorput(int _base = 2, myLong _currentNumber = 1);
	virtual double Generate();


private:
	int base;

};

std::vector<int> IntToInversePAdicExpansion(int n, int base);