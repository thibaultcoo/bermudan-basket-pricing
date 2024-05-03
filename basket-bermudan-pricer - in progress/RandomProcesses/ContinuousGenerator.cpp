#include "pch.h"
#include "ContinuousGenerator.h"


ContinuousGenerator::ContinuousGenerator()
{
}

ContinuousGenerator::ContinuousGenerator(UniformGenerator* _gen)
	:generator(_gen)
{
}