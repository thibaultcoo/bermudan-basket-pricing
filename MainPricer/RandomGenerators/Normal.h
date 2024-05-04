#pragma once
#include "ContinuousGenerator.h"

class Normal : public ContinuousGenerator
{
	protected:
		double mean;
		double var;

	public:
		Normal();
		Normal(double _mean, double _var, UniformGenerator* _gen);
};

class NormalBoxMuller : public Normal
{
	protected:
		bool NewSimulation = true;
		double SecondNormal;

	public:
		NormalBoxMuller();
		NormalBoxMuller(double _mean, double _var, UniformGenerator* _gen);
		double Generate();
};

class NormalCLT : public Normal
{
	public:
		NormalCLT();
		NormalCLT(double _mean, double _var, UniformGenerator* _gen);
		double Generate();
};

class NormalInvCdf :public Normal
{
	public:
		NormalInvCdf();
		NormalInvCdf(double _mean, double _var, UniformGenerator* _gen);
		virtual double Generate();
};

// Helper Functions
double RationalApproximation(double t);
double NormalCDFInverse(double p);