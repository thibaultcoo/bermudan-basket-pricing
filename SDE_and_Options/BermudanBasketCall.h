#pragma once
#include "Option.h"
#include <eigen-3.4.0/Eigen/Dense>

//Bermudan basket Call pricing using the Longstaff-Shawratz algorithm

class Bermudan : public Option
{
public:
	Bermudan();
	Bermudan(StochasticProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights, std::vector<double> _exeDates, int _L = 5);
	Bermudan(StochasticProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights, std::vector<double> _exeDates, std::vector<double> _S, Eigen::MatrixXd _VarCovar, int _L = 5);
	double ComputePrice(int NbSim, bool antithetic = false);
	double ComputePrice_ControlVariate(int NbSim);
	double ComputePrice_VDC(int NbSim);

private:
	std::vector<double> exeDates;
	int L;
	std::vector<double> weights;
	std::vector<double> S;
	Eigen::MatrixXd VarCovar;
};

double Compute_E_Ybis(std::vector<double> S, std::vector<double> weights, double K, double r, Eigen::MatrixXd VarCovar, double T);

