#include "pch.h"
#include "Normal.h"


Normal::Normal()
{
}

Normal::Normal(double _mean, double _var, UniformGenerator* _gen)
	:mean(_mean), var(_var), ContinuousGenerator(_gen)
{
}