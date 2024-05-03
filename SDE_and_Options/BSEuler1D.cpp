#include "pch.h"
#include "BSEuler1D.h"
#include "SinglePath.h"
#include "iostream"
#include "VanDerCorput.h"
#include "NormalInvCdf.h"
#include <algorithm>
#include <random>

BSEuler1D::BSEuler1D()
{
}

BSEuler1D::BSEuler1D(RandomGenerator* _gen, double _s, double _r, double _vol)
	: BlackScholes1D(_gen, _s, _r, _vol)
{
}

void BSEuler1D::Simulate(double start_time, double end_time, size_t nb_steps)
{
	SinglePath* path = new SinglePath(start_time, end_time, nb_steps);
	path->AddValue(s);
	double last = s;
	double dt = (end_time - start_time) / nb_steps;

	for (size_t i = 0; i < nb_steps; i++)
	{
		double dW = pow(dt, 0.5) * gen->Generate();
		double next = last + last * (r * dt + vol * dW);
		path->AddValue(next);
		last = next;
	}

	if (paths[0])
	{
		delete paths[0];
	}
	paths[0] = path;

}

void BSEuler1D::Simulate_Antithetic(double start_time, double end_time, size_t nb_steps)
{
	SinglePath* path = new SinglePath(start_time, end_time, nb_steps);
	SinglePath* path_anti = new SinglePath(start_time, end_time, nb_steps);
	path->AddValue(s);
	path_anti->AddValue(s);
	double last = s;
	double last_anti = s;
	double dt = (end_time - start_time) / nb_steps;
	for (size_t i = 0; i < nb_steps; i++)
	{
		double dW = pow(dt, 0.5) * gen->Generate();
		double next = last + last * (r * dt + vol * dW);
		double next_anti = last_anti + last_anti * (r * dt - vol * dW);
		path->AddValue(next);
		path_anti->AddValue(next_anti);
		last = next;
		last_anti = next_anti;
	}

	if (paths[0])
	{
		delete paths[0];
	}
	paths[0] = path;

	if (paths_antithetic[0])
	{
		delete paths_antithetic[0];
	}
	paths_antithetic[0] = path_anti;

}

void BSEuler1D::Simulate_VDC(double start_time, double end_time, size_t nb_steps, myLong sim, myLong nbSim)
{
	SinglePath* path = new SinglePath(start_time, end_time, nb_steps);
	path->AddValue(s);
	double last = s;
	double dt = (end_time - start_time) / nb_steps;
	std::vector<double> x;

	if (v_VDC.size() != nbSim)
	{
		for (size_t i = 0; i < nbSim; i++)
		{
			v_VDC.push_back(gen->Generate());
		}

	}
	double dW = pow(dt, 0.5) * v_VDC[sim];
	double next = last + last * (r * dt + vol * dW);
	path->AddValue(next);
	last = next;

	for (unsigned int i = 1; i < nb_steps; i++)
	{
		x = v_VDC;
		auto rng = std::default_random_engine{ i };
		std::shuffle(x.begin(), x.end(), rng); //Needed permutation for quasi-monte carlo
		double dW = pow(dt, 0.5) * x[sim];
		double next = last + last * (r * dt + vol * dW);
		path->AddValue(next);
		last = next;
	}
	if (paths[0])
	{
		delete paths[0];
	}
	paths[0] = path;
}

