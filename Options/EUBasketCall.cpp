#include "EUBasketCall.h"
#include <numeric>

EUBasketCall::EUBasketCall() { }

EUBasketCall::EUBasketCall(StochasticProcess* process, double strike, std::vector<double> rates, double maturity, std::vector<double> weights)
	: Options(process, strike, rates, maturity), weights(weights) { }

EUBasketCall::EUBasketCall(StochasticProcess* process, double strike, std::vector<double> rates, double maturity, std::vector<double> weights, std::vector<double> spots, Eigen::MatrixXd corrMatrix)
	: Options(process, strike, rates, maturity), weights(weights), spots(spots), corrMatrix(corrMatrix) { }

double EUBasketCall::price(int nbSim)
{
	double sum = 0.;
	double path_price, weighted_spots;

	paths_prices.clear();
	paths_prices.resize(nbSim);

	for (int n = 0; n < nbSim; ++n)
	{
		process->Simulate(0, maturity, maturity * 365);
		std::vector<double> final_spots = process->Get_ValueND(maturity);
		weighted_spots = std::inner_product(std::begin(weights), std::end(weights), std::begin(final_spots), 0.0);
		path_price = std::max(weighted_spots - strike, 0.);
		sum += path_price;
		paths_prices[n] = path_price;
	}

	double price = std::exp(-rates[0] * maturity) * (sum / nbSim);

	return price;
}

double EUBasketCall::priceAntithetic(int nbSim)
{
	double sum = 0.;
	double path_price, weighted_spots;
	double path_price_anti, weighted_spots_anti;

	paths_prices.clear();
	paths_prices.resize(nbSim);

	for (int n = 0; n < nbSim; ++n)
	{
		process->Simulate_Antithetic(0, maturity, maturity * 365);
		std::vector<double> final_spots = process->Get_ValueND(maturity);
		std::vector<double> final_spots_anti = process->Get_ValueND_antithetic(maturity);
		weighted_spots = std::inner_product(std::begin(weights), std::end(weights), std::begin(final_spots), 0.0);
		weighted_spots_anti = std::inner_product(std::begin(weights), std::end(weights), std::begin(final_spots_anti), 0.0);
		path_price = std::max(weighted_spots - strike, 0.);
		path_price_anti = std::max(weighted_spots_anti - strike, 0.);
		sum += (path_price + path_price_anti);
		paths_prices[n] = (path_price + path_price_anti) / 2;
	}

	double price = std::exp(-rates[0] * maturity) * (sum / (nbSim * 2));

	return price;
}

double EUBasketCall::priceControlVariate(int nbSim)
{
	double sum = 0.;
	double weighted_spots, weighted_log_spots;
	double payoff, control_adjusted_payoff, log_payoff;

	paths_prices.clear();
	paths_prices.resize(nbSim);

	double control_variate_payoff = compute_expected_value_control_variate(spots, weights, strike, rates[0], corrMatrix, maturity);

	for (int n = 0; n < nbSim; ++n)
	{
		process->Simulate(0, maturity, maturity * 365);
		std::vector<double> final_spots = process->Get_ValueND(maturity);
		std::vector<double> log_spots(spots.size());
		for (int j = 0; j < spots.size(); j++)
			log_spots[j] = log(final_spots[j]);

		weighted_spots = std::inner_product(std::begin(weights), std::end(weights), std::begin(final_spots), 0.0);
		weighted_log_spots = std::inner_product(std::begin(weights), std::end(weights), std::begin(log_spots), 0.0);

		payoff = std::max(weighted_spots - strike, 0.);
		log_payoff = std::max(exp(weighted_log_spots) - strike, 0.);
		control_adjusted_payoff = payoff - log_payoff + control_variate_payoff;

		sum += control_adjusted_payoff;
		paths_prices[n] = control_adjusted_payoff;
	}

	double price = std::exp(-rates[0] * maturity) * (sum / nbSim);

	return price;
}

double EUBasketCall::priceVDC(int nbSim)
{
	double sum = 0.;
	double path_price, weighted_spots, weighted_spots_anti;

	paths_prices.clear();
	paths_prices.resize(nbSim);

	for (int n = 0; n < nbSim; ++n)
	{
		process->Simulate_VDC(0, maturity, maturity * 365, n, nbSim);
		std::vector<double> final_spots = process->Get_ValueND(maturity);
		weighted_spots = std::inner_product(std::begin(weights), std::end(weights), std::begin(final_spots), 0.0);
		path_price = std::max(weighted_spots - strike, 0.);
		sum += path_price;
		paths_prices[n] = path_price;
	}
	
	double price = std::exp(-rates[0] * maturity) * (sum / nbSim);

	return price;
}

// Functions
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

