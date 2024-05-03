#include "Tools.h"
#include <cmath>
#include <stdexcept>

// Compute the cumulative distribution function for a standard normal distribution
double norm_cdf(double x)
{
	return 0.5 * std::erfc(-x / std::sqrt(2));
}

// Function to calculate a call option price using Black-Scholes closed formula
double BS_Call(double spot, double strike, double volatility, double maturity, double rate)
{
	double sigma = volatility * std::sqrt(maturity);
	double d1 = (std::log(spot / strike) + (rate + 0.5 * pow(volatility, 2)) * maturity) / sigma;
	double d2 = d1 - sigma;

	return spot * norm_cdf(d1) - strike * std::exp(-rate * maturity) * norm_cdf(d2);;
}

// Function to compute the expected value using the control variate method
double compute_expected_value_control_variate(std::vector<double> spots, std::vector<double> weights, double strike, double rate, Eigen::MatrixXd corrMatrix, double maturity) 
{
	Eigen::VectorXd weightsVec = Eigen::Map<Eigen::VectorXd>(weights.data(), weights.size());
	Eigen::MatrixXd choleskyDecomp;

	// Check if the matrix is singular
	if (corrMatrix.determinant() == 0) 
	{
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(corrMatrix);
		choleskyDecomp = es.eigenvectors() * es.eigenvalues().cwiseMax(0).cwiseSqrt().asDiagonal();
	}
	else 
	{
		// Cholesky decomposition for non-singular matrices
		choleskyDecomp = corrMatrix.llt().matrixL();
	}

	double spotProduct = 1.0;
	for (size_t i = 0; i < spots.size(); ++i)
		spotProduct *= std::pow(spots[i], weights[i]);

	Eigen::VectorXd variances = (choleskyDecomp * choleskyDecomp.transpose()).colwise().sum();
	double adjustedRate = rate - 0.5 * weightsVec.dot(variances) + 0.5 * weightsVec.transpose() * choleskyDecomp * choleskyDecomp.transpose() * weightsVec;
	double volAdjusted = std::sqrt(weightsVec.transpose() * choleskyDecomp * choleskyDecomp.transpose() * weightsVec);

	return BS_Call(spotProduct, strike, volAdjusted, maturity, adjustedRate);
}