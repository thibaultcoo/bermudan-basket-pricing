#include "pch.h"
#include "EuropeanCall.h"
#include "BSEuler1D.h"
#include "BSMilstein1D.h"
#include "SinglePath.h"
#include <iostream>     
#include <algorithm> 
#include <cmath>

EuropeanCall::EuropeanCall()
{
}

EuropeanCall::EuropeanCall(StochasticProcess* _process, double _K, std::vector<double> _r, double _T)
	: Option(_process, _K, _r, _T)
{

}

double EuropeanCall::ComputePrice(int NbSim, bool antithetic)
{
	double somme = 0.;
	double last_value1;
	double last_value2;
	double price;

	v.clear();
	v.resize(NbSim);

	if (antithetic)
	{
		for (int n = 0; n < NbSim; ++n)
		{
			process->Simulate_Antithetic(0, T, T * 365);
			last_value1 = std::max(process->Get_Value(T) - K, 0.);
			last_value2 = std::max(process->Get_Value_antithetic(T) - K, 0.);
			somme = somme + last_value1 + last_value2;
			v[n] = (last_value1 + last_value2) / 2;
		}
		price = std::exp(-r[0] * T) * (somme / (NbSim * 2));
	}
	else
	{
		for (int n = 0; n < NbSim; ++n)
		{
			process->Simulate(0, T, T * 365);
			last_value1 = std::max(process->Get_Value(T) - K, 0.);
			somme = somme + last_value1;
			v[n] = last_value1;
		}
		price = std::exp(-r[0] * T) * (somme / NbSim);
	}

	return price;

}

double EuropeanCall::ComputePrice_ControlVariate(int NbSim)
{
	double somme = 0.;
	double last_value;
	double price;

	v.clear();
	v.resize(NbSim);

	for (int n = 0; n < NbSim; ++n)
	{
		process->Simulate(0, T, T * 365);
		last_value = std::max(K - process->Get_Value(T), 0.) + std::exp(r[0] * T) * process->Get_Value(0) - K;
		somme = somme + last_value;
		v[n] = last_value;
	}

	price = std::exp(-r[0] * T) * (somme / NbSim);

	return price;
}

double EuropeanCall::ComputePrice_VDC(int NbSim)
{
	double somme = 0.;
	double last_value1;
	double last_value2;
	double price;

	v.clear();
	v.resize(NbSim);

	for (int n = 0; n < NbSim; ++n)
	{
		process->Simulate_VDC(0, T, T * 365, n, NbSim);
		last_value1 = std::max(process->Get_Value(T) - K, 0.);
		somme = somme + last_value1;
		v[n] = last_value1;
	}
	price = std::exp(-r[0] * T) * (somme / NbSim);
	return price;
}