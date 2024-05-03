#include "pch.h"
#include "FiniteSet.h"

FiniteSet::FiniteSet()
{
}

FiniteSet::FiniteSet(UniformGenerator* _gen, std::vector<double> _values, std::vector<double> _probas)
	:DiscreteGenerator(_gen), values(_values), probas(_probas)
{
}

double FiniteSet::Generate()
{
	std::vector<double> SumProbas(probas.size());
	SumProbas[0] = probas[0];
	for (int i = 1; i < probas.size(); i++)
	{
		SumProbas[i] = SumProbas[i - 1] + probas[i];
	}
	double u = generator->Generate();

	for (int i = 0; i < probas.size(); i++)
	{
		if (u < SumProbas[i])
		{
			return values[i];
		}
	}

}