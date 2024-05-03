#include "pch.h"
#include "BlackScholesnD.h"

BlackScholesnD::BlackScholesnD()
{
}

BlackScholesnD::BlackScholesnD(RandomGenerator* _gen, std::vector<double> _s, std::vector<double> _r, Eigen::MatrixXd _VarCovar, int dim)
	: StochasticProcess(_gen, dim), s(_s), r(_r), VarCovar(_VarCovar)
{
}


