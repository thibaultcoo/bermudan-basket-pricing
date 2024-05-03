#include "pch.h"
#include "StochasticProcess.h"
#include <vector>
#include "SinglePath.h"

StochasticProcess::StochasticProcess()
{
}

StochasticProcess::StochasticProcess(RandomGenerator* _gen, int _dim)
	:gen(_gen), dim(_dim), paths(_dim, new SinglePath()), paths_antithetic(_dim, new SinglePath())
{
}

void StochasticProcess::add_path(SinglePath* Path)
{
	paths.push_back(Path);
}

SinglePath* StochasticProcess::GetPath(int dimension)
{
	return paths[0];
}

const double StochasticProcess::Get_Value(double time, int dim)
{
	return paths[dim]->GetValue(time);
}

const double StochasticProcess::Get_Value_antithetic(double time, int dim)
{
	return paths_antithetic[dim]->GetValue(time);
}


const std::vector<double> StochasticProcess::Get_ValueND(double time)
{
	int d = dim;

	std::vector<double> res(d);
	for (int i = 0; i < d; i++)
	{
		res[i] = Get_Value(time, i);
	}
	return res;
}

const std::vector<double> StochasticProcess::Get_ValueND_antithetic(double time)
{
	int d = dim;

	std::vector<double> res(d);
	for (int i = 0; i < d; i++)
	{
		res[i] = Get_Value_antithetic(time, i);
	}
	return res;
}


