#include "pch.h"
#include "RandomGenerator.h"
#include <math.h>
#include "VanDerCorput.h"
#include "NormalInvCdf.h"
#include <iostream>

RandomGenerator::RandomGenerator()
{
}

//Get the mean of any random processes implemented
double RandomGenerator::Mean(myLong nbSim)
{
	double result = 0.;
	lastGeneratedNumbers = std::vector<double>(nbSim);

	for (myLong i = 0; i < nbSim; i++)
	{
		double currentNumber = Generate();
		result += currentNumber / nbSim;
		lastGeneratedNumbers[i] = currentNumber;
	}

	return result;
}

//Get the variance of any random processes implemented
double RandomGenerator::Variance(myLong nbSim)
{
	double result = 0.;
	double mean = Mean(nbSim);

	for (myLong i = 0; i < nbSim; i++)
	{
		double currentNumber = lastGeneratedNumbers[i];
		result += pow(currentNumber - mean, 2) / nbSim;
	}

	return result;
}

//Generate vector with random values
std::vector<double> RandomGenerator::GenerateVector(int dimension)
{
	std::vector<double> result(dimension);
	for (int i = 0; i < dimension; i++)
	{
		result[i] = Generate();
	}
	return result;
}

//Generate vector with random values using Van Der Corupt Sequence
std::vector<double> RandomGenerator::GenerateVectorVDC(int dimension, myLong sim)
{
	std::vector<double> result(dimension);
	std::vector<double> prime_numbers = first_primeNumbers(dimension);
	result[0] = Generate();
	for (int i = 1; i < dimension; i++)
	{
		double prime_number = prime_numbers[i];
		VanDerCorput* Vdc = new VanDerCorput(prime_number, sim + 1);
		//std::cout << "prime_number : " << prime_number << std::endl;
		NormalInvCdf* Normal_VDC = new NormalInvCdf(0, 1, Vdc);
		result[i] = Normal_VDC->Generate();
	}
	return result;
}

//Get the N first prime numbers
std::vector<double> first_primeNumbers(myLong N)
{
	std::vector<double> results;
	myLong counter;
	myLong num;
	myLong i = 3;
	results.push_back(2);

	for (counter = 2; counter <= N; i++)
	{
		for (num = 2; num < i; num++)
		{
			if (i % num == 0)
				break;
		}
		if (num == i)
		{
			results.push_back(i);
			counter++;
		}
	}

	return results;
}