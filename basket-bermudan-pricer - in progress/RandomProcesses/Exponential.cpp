#include "pch.h"
#include "Exponential.h"

Exponential::Exponential()
{
}

Exponential::Exponential(double _lambda, UniformGenerator* _gen)
	:ContinuousGenerator(_gen), lambda(_lambda)
{
}