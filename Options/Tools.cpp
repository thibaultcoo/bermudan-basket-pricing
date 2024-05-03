#include "Tools.h"

double norm_cdf(double x)
{
	return 0.5 * std::erfc(-x / std::sqrt(2));
}

double BS_Call(double spot, double strike, double volatility, double maturity, double rate)
{
	double sigma = volatility * std::sqrt(maturity);
	double d1 = (std::log(spot / strike) + (rate + 0.5 * pow(volatility, 2)) * maturity) / sigma;
	double d2 = d1 - sigma;

	double price = spot * norm_cdf(d1) - strike * std::exp(-rate * maturity) * norm_cdf(d2);

	return price;
}

double compute_expected_value_control_variate(std::vector<double> S, std::vector<double> weights, double K, double r, Eigen::MatrixXd VarCovar, double T)
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

	double price = BS_Call(S_Y, K, vol_Y, T, R_Y);

	return price;
}