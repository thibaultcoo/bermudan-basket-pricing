#include "pch.h"
#include "NormalCLT.h"

NormalCLT::NormalCLT()
{
}

NormalCLT::NormalCLT(double _mean, double _var, UniformGenerator* _gen)
	:Normal(_mean, _var, _gen)
{
}

double NormalCLT::Generate()
{
	double sumUniform = 0;
	for (int i = 0; i < 12; i++)
	{
		sumUniform += generator->Generate();
	}

	return (sumUniform - 6) * pow(var, 0.5) + mean;
}