#pragma once
#include "QuasiGenerator.h"

class VanDerCorput : public QuasiGenerator
{
	protected:
		int base;

	public:
		VanDerCorput(int _base = 2, myLong _currentNumber = 1);
		virtual double Generate();
};

// Helper function
std::vector<int> InversePAdicExpansion(int n, int base);