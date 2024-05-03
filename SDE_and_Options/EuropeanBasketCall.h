#pragma once
#include "Option.h"
#include <eigen-3.4.0/Eigen/Dense>

//European basket call pricing

class European :public Option
{
public:
	European();
	European(StochasticProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights);
	European(StochasticProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights, std::vector<double> _S, Eigen::MatrixXd _VarCovar);
	double ComputePrice(int NbSim, bool antithetic = false);
	double ComputePrice_ControlVariate(int NbSim);
	double ComputePrice_VDC(int NbSim);

private:
	std::vector<double> weights;
	std::vector<double> S;
	Eigen::MatrixXd VarCovar;
};

double Compute_E_Y(std::vector<double> S, std::vector<double> weights, double K, double r, Eigen::MatrixXd VarCovar, double T);