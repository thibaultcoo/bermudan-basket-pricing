#include "pch.h"
#include "EuropeanBasketCall.h"
#include <numeric>
#include <iterator>
#include "ClosedFormulas.h"

European::European()
{
}

European::European(StochasticProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights)
	: Option(_process, _K, _r, _T), weights(_weights)
{

}
European::European(StochasticProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights, std::vector<double> _S, Eigen::MatrixXd _VarCovar)
	: Option(_process, _K, _r, _T), weights(_weights), S(_S), VarCovar(_VarCovar)
{
}

double European::ComputePrice(int NbSim, bool antithetic)
{
	double somme = 0.;
	double last_value1;
	double last_value2;
	double WS_T;
	double WS_T_anti;
	double price;

	v.clear();
	v.resize(NbSim);

	if (antithetic)
	{
		for (int n = 0; n < NbSim; ++n)
		{
			process->Simulate_Antithetic(0, T, T * 365);
			std::vector<double> S_T = process->Get_ValueND(T);
			std::vector<double> S_T_anti = process->Get_ValueND_antithetic(T);
			WS_T = std::inner_product(std::begin(weights), std::end(weights), std::begin(S_T), 0.0);
			WS_T_anti = std::inner_product(std::begin(weights), std::end(weights), std::begin(S_T_anti), 0.0);
			last_value1 = std::max(WS_T - K, 0.);
			last_value2 = std::max(WS_T_anti - K, 0.);
			somme = somme + last_value1 + last_value2;
			v[n] = (last_value1 + last_value2) / 2;
		}
		price = std::exp(-r[0] * T) * (somme / (NbSim * 2));
	}
	else
	{
		for (int n = 0; n < NbSim; ++n)
		{
			process->Simulate(0, T, T * 365);
			std::vector<double> S_T = process->Get_ValueND(T);
			WS_T = std::inner_product(std::begin(weights), std::end(weights), std::begin(S_T), 0.0);
			last_value1 = std::max(WS_T - K, 0.);
			somme = somme + last_value1;
			v[n] = last_value1;
		}
		price = std::exp(-r[0] * T) * (somme / NbSim);
	}

	return price;

}

double European::ComputePrice_ControlVariate(int NbSim)
{
	double somme = 0.;
	double last_value;
	double WS_T;
	double WS_T_L;
	double X;
	double Y;
	double X_prime;

	v.clear();
	v.resize(NbSim);

	double E_Y = Compute_E_Y(S, weights, K, r[0], VarCovar, T);

	for (int n = 0; n < NbSim; ++n)
	{
		process->Simulate(0, T, T * 365);
		std::vector<double> S_T = process->Get_ValueND(T);
		std::vector<double> S_T_L(S.size());
		for (int j = 0; j < S.size(); j++)
		{
			S_T_L[j] = log(S_T[j]);
		}
		WS_T = std::inner_product(std::begin(weights), std::end(weights), std::begin(S_T), 0.0);
		WS_T_L = std::inner_product(std::begin(weights), std::end(weights), std::begin(S_T_L), 0.0);

		X = std::max(WS_T - K, 0.);
		Y = std::max(exp(WS_T_L) - K, 0.);
		X_prime = X - Y + E_Y;

		somme = somme + X_prime;
		v[n] = X_prime;
	}

	double price = std::exp(-r[0] * T) * (somme / NbSim);
	return price;
}

double Compute_E_Y(std::vector<double> S, std::vector<double> weights, double K, double r, Eigen::MatrixXd VarCovar, double T)
{

	double S_Y = 1;
	double R_Y;
	double vol_Y;

	Eigen::VectorXd weights_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());
	Eigen::MatrixXd B;

	Eigen::EigenSolver<Eigen::MatrixXd> es(VarCovar);

	if (VarCovar.determinant() == 0)
	{
		Eigen::MatrixXd D = es.pseudoEigenvalueMatrix();
		Eigen::MatrixXd P = es.pseudoEigenvectors();
		B = P * D;
	}
	else
	{
		B = VarCovar.llt().matrixL();
	}

	for (int i = 0; i < S.size(); i++)
	{
		S_Y *= pow(S[i], weights[i]);
	}

	Eigen::MatrixXd sigma2 = B * B;
	Eigen::VectorXd sigmai = sigma2.colwise().sum();
	double R1 = 0.5 * weights_M.transpose() * sigmai;
	double R2 = 0.5 * weights_M.transpose() * B * B.transpose() * weights_M;
	R_Y = r - R1 + R2;

	vol_Y = pow(weights_M.transpose() * B * B.transpose() * weights_M, 0.5);

	double price = bs_price_call(S_Y, K, vol_Y, T, R_Y);

	return price;

}

double European::ComputePrice_VDC(int NbSim)
{
	double somme = 0.;
	double last_value1;
	double last_value2;
	double WS_T;
	double WS_T_anti;
	double price;

	v.clear();
	v.resize(NbSim);

	for (int n = 0; n < NbSim; ++n)
	{
		process->Simulate_VDC(0, T, T * 365, n, NbSim);
		std::vector<double> S_T = process->Get_ValueND(T);
		WS_T = std::inner_product(std::begin(weights), std::end(weights), std::begin(S_T), 0.0);
		last_value1 = std::max(WS_T - K, 0.);
		somme = somme + last_value1;
		v[n] = last_value1;
	}
	price = std::exp(-r[0] * T) * (somme / NbSim);

	return price;;
}