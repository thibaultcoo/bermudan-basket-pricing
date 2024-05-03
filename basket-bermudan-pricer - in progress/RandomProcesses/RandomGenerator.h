#pragma once
#include <vector>

typedef unsigned long long myLong;

class RandomGenerator
{
	public:
		// Constructor
		RandomGenerator();

		// Methods
		virtual double Generate() = 0;
		std::vector<double> GenerateVector(int dimension);
		std::vector<double> GenerateVectorVDC(int dimension, myLong sim);
		double Mean(myLong nbSim);
		double Variance(myLong nbSim);

		// Members
		std::vector <double> lastGeneratedNumbers;
};

std::vector<double> first_primeNumbers(myLong N);