#include "pch.h"
#include "BlackScholes1D.h"

BlackScholes1D::BlackScholes1D()
{
}

BlackScholes1D::BlackScholes1D(RandomGenerator* _gen, double _s, double _r, double _vol)
	: StochasticProcess(_gen, 1), s(_s), r(_r), vol(_vol)
{
}