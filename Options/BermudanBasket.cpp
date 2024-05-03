#include "BermudanBasket.h"
#include "Tools.h"
#include <algorithm>
#include <cmath>

// Default Constructor
BermudanBasket::BermudanBasket() { }

// Initializer Constructors
BermudanBasket::BermudanBasket(StochasticProcess* process, double strike, std::vector<double> rates, double maturity, std::vector<double> weights, std::vector<double> exerciseDates, int L)
	: Options(process, strike, rates, maturity), exerciseDates(exerciseDates), L(L), weights(weights) { }

BermudanBasket::BermudanBasket(StochasticProcess* process, double strike, std::vector<double> rates, double maturity, std::vector<double> weights, std::vector<double> exerciseDates, std::vector<double> spots, Eigen::MatrixXd corrMatrix, int L)
	: Options(process, strike, rates, maturity), exerciseDates(exerciseDates), L(L), weights(weights), spots(spots), corrMatrix(corrMatrix) { }

// Pricing Methods
double BermudanBasket::price(int nbSim)
{
	int numExerciseDates = exerciseDates.size();
	Eigen::MatrixXd basketValues(nbSim, numExerciseDates);
	Eigen::MatrixXd exerciseTimes(nbSim, numExerciseDates);
	Eigen::VectorXd basketWeightsVector = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

	double currentExerciseTime;
	double optionPrice = 0.0;

	paths_prices.clear();
	paths_prices.resize(nbSim);

	// Simulate paths and calculate basket values at each exercise date
	for (int simIndex = 0; simIndex < nbSim; ++simIndex)
	{
		exerciseTimes(simIndex, numExerciseDates - 1) = maturity;
		process->Simulate(0, maturity, maturity * 365);

		for (int dateIndex = 0; dateIndex < numExerciseDates; dateIndex++)
		{
			std::vector<double> spotPrices = process->Get_ValueND(exerciseDates[dateIndex]);
			Eigen::VectorXd spotPricesVector = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(spotPrices.data(), spotPrices.size());
			basketValues(simIndex, dateIndex) = basketWeightsVector.transpose() * spotPricesVector;
		}
	}

	// Backward induction loop
	for (int dateIndex = numExerciseDates - 2; dateIndex >= 0; dateIndex--)
	{
		backwardInduction(dateIndex, nbSim, basketValues, exerciseTimes);
	}

	// Calculate option price
	for (int simIndex = 0; simIndex < nbSim; simIndex++)
	{
		optionPrice += discountedPayoff(simIndex, numExerciseDates, basketValues, exerciseTimes);
	}

	optionPrice /= nbSim;

	return optionPrice;
}

double BermudanBasket::priceAntithetic(int numSimulations)
{
	int numExerciseDates = exerciseDates.size();
	Eigen::MatrixXd basketValuesRegular(numSimulations, numExerciseDates);
	Eigen::MatrixXd basketValuesAntithetic(numSimulations, numExerciseDates);
	Eigen::MatrixXd exerciseTimesRegular(numSimulations, numExerciseDates);
	Eigen::MatrixXd exerciseTimesAntithetic(numSimulations, numExerciseDates);
	Eigen::VectorXd basketWeightsVector = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

	double sumRegular = 0.0;
	double sumAntithetic = 0.0;

	paths_prices.clear();
	paths_prices.resize(numSimulations);

	// Simulate paths and calculate basket values at each exercise date for regular and antithetic paths
	for (int simIndex = 0; simIndex < numSimulations; ++simIndex)
	{
		exerciseTimesRegular(simIndex, numExerciseDates - 1) = maturity;
		exerciseTimesAntithetic(simIndex, numExerciseDates - 1) = maturity;
		process->Simulate_Antithetic(0, maturity, maturity * 365);

		for (int dateIndex = 0; dateIndex < numExerciseDates; dateIndex++)
		{
			std::vector<double> spotPricesRegular = process->Get_ValueND(exerciseDates[dateIndex]);
			Eigen::VectorXd spotPricesVectorRegular = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(spotPricesRegular.data(), spotPricesRegular.size());
			basketValuesRegular(simIndex, dateIndex) = basketWeightsVector.transpose() * spotPricesVectorRegular;

			std::vector<double> spotPricesAntithetic = process->Get_ValueND_antithetic(exerciseDates[dateIndex]);
			Eigen::VectorXd spotPricesVectorAntithetic = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(spotPricesAntithetic.data(), spotPricesAntithetic.size());
			basketValuesAntithetic(simIndex, dateIndex) = basketWeightsVector.transpose() * spotPricesVectorAntithetic;
		}
	}

	// Perform backward induction for regular and antithetic paths
	for (int dateIndex = numExerciseDates - 2; dateIndex >= 0; dateIndex--)
	{
		backwardInduction(dateIndex, numSimulations, basketValuesRegular, exerciseTimesRegular);
		backwardInduction(dateIndex, numSimulations, basketValuesAntithetic, exerciseTimesAntithetic);
	}

	// Calculate option price for regular and antithetic paths
	for (int simIndex = 0; simIndex < numSimulations; simIndex++)
	{
		sumRegular += discountedPayoff(simIndex, numExerciseDates, basketValuesRegular, exerciseTimesRegular);
		sumAntithetic += discountedPayoff(simIndex, numExerciseDates, basketValuesAntithetic, exerciseTimesAntithetic);
	}

	double price = ((sumRegular + sumAntithetic) / (2 * numSimulations));

	return price;
}

double BermudanBasket::priceControlVariate(int numSimulations)
{
	int numExerciseDates = exerciseDates.size();
	Eigen::MatrixXd basketValues(numSimulations, numExerciseDates);
	Eigen::MatrixXd basketValuesLog(numSimulations, numExerciseDates);
	Eigen::MatrixXd exerciseTimes(numSimulations, numExerciseDates);
	Eigen::VectorXd basketWeightsVector = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

	double sum = 0.0;
	double E_Y = compute_expected_value_control_variate(spots, weights, strike, rates[0], corrMatrix, maturity);

	paths_prices.clear();
	paths_prices.resize(numSimulations);

	// Simulate paths and calculate basket values at each exercise date
	for (int simIndex = 0; simIndex < numSimulations; ++simIndex)
	{
		exerciseTimes(simIndex, numExerciseDates - 1) = maturity;
		process->Simulate(0, maturity, maturity * 365);

		for (int dateIndex = 0; dateIndex < numExerciseDates; dateIndex++)
		{
			std::vector<double> spotPrices = process->Get_ValueND(exerciseDates[dateIndex]);
			std::vector<double> logSpotPrices(spotPrices.size());
			std::transform(spotPrices.begin(), spotPrices.end(), logSpotPrices.begin(), [](double price) { return log(price); });

			Eigen::VectorXd spotPricesVector = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(spotPrices.data(), spotPrices.size());
			Eigen::VectorXd logSpotPricesVector = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(logSpotPrices.data(), logSpotPrices.size());

			basketValues(simIndex, dateIndex) = basketWeightsVector.transpose() * spotPricesVector;
			basketValuesLog(simIndex, dateIndex) = basketWeightsVector.transpose() * logSpotPricesVector;
		}
	}

	// Perform backward induction
	for (int dateIndex = numExerciseDates - 2; dateIndex >= 0; dateIndex--)
	{
		double currentExerciseTime = exerciseDates[dateIndex];
		Eigen::MatrixXd regressionMatrix(numSimulations, L);
		Eigen::MatrixXd regressionTarget(numSimulations, 1);

		// Calculate regression coefficients
		for (int simIndex = 0; simIndex < numSimulations; simIndex++)
		{
			double nextExerciseTime = exerciseTimes(simIndex, dateIndex + 1);
			int lastIndex = numExerciseDates - 1;
			while (nextExerciseTime != exerciseDates[lastIndex])
				lastIndex--;

			double payoffDiff = std::max(basketValues(simIndex, lastIndex) - strike, 0.) - std::max(exp(basketValuesLog(simIndex, lastIndex)) - strike, 0.) + E_Y;

			regressionTarget(simIndex, 0) = std::exp(-rates[0] * (exerciseTimes(simIndex, dateIndex + 1) - currentExerciseTime)) * payoffDiff;

			for (int powerIndex = 0; powerIndex < L; powerIndex++)
				regressionMatrix(simIndex, powerIndex) = pow(basketValues(simIndex, dateIndex), powerIndex);
		}

		// Solve for regression coefficients
		Eigen::VectorXd regressionCoefficients(regressionMatrix.colPivHouseholderQr().solve(regressionTarget));

		// Determine optimal exercise strategy
		for (int simIndex = 0; simIndex < numSimulations; simIndex++)
		{
			double actualPayoffDiff = std::max(basketValues(simIndex, dateIndex) - strike, 0.) - std::max(exp(basketValuesLog(simIndex, dateIndex)) - strike, 0.) + E_Y;
			double expectedPayoffDiff = regressionMatrix.row(simIndex) * regressionCoefficients;

			exerciseTimes(simIndex, dateIndex) = (actualPayoffDiff >= expectedPayoffDiff) ? currentExerciseTime : exerciseTimes(simIndex, dateIndex + 1);
		}
	}

	// Calculate option price
	for (int simIndex = 0; simIndex < numSimulations; simIndex++)
	{
		// Find exercise time at initial date
		double initialExerciseTime = exerciseTimes(simIndex, 0);
		int lastIndex = numExerciseDates - 1;

		while (initialExerciseTime != exerciseDates[lastIndex])
			lastIndex--;

		double discountedPayoff = std::exp(-rates[0] * initialExerciseTime) * (std::max(basketValues(simIndex, lastIndex) - strike, 0.) - std::max(exp(basketValuesLog(simIndex, lastIndex)) - strike, 0.) + E_Y);

		paths_prices[simIndex] = discountedPayoff;
		sum += discountedPayoff;
	}

	double price = (sum / numSimulations);

	return price;
}

double BermudanBasket::priceVDC(int numSimulations)
{
	int numExerciseDates = exerciseDates.size();
	Eigen::MatrixXd basketValues(numSimulations, numExerciseDates);
	Eigen::MatrixXd exerciseTimes(numSimulations, numExerciseDates);
	Eigen::VectorXd basketWeightsVector = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

	double sum = 0.0;

	paths_prices.clear();
	paths_prices.resize(numSimulations);

	// Simulate paths and calculate basket values at each exercise date
	for (int simIndex = 0; simIndex < numSimulations; ++simIndex)
	{
		exerciseTimes(simIndex, numExerciseDates - 1) = maturity;
		process->Simulate_VDC(0, maturity, maturity * 365, simIndex, numSimulations);

		for (int dateIndex = 0; dateIndex < numExerciseDates; dateIndex++)
		{
			std::vector<double> spotPrices = process->Get_ValueND(exerciseDates[dateIndex]);
			Eigen::VectorXd spotPricesVector = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(spotPrices.data(), spotPrices.size());
			basketValues(simIndex, dateIndex) = basketWeightsVector.transpose() * spotPricesVector;
		}
	}

	// Perform backward induction
	for (int dateIndex = numExerciseDates - 2; dateIndex >= 0; dateIndex--)
	{
		backwardInduction(dateIndex, numSimulations, basketValues, exerciseTimes);
	}

	// Calculate option price
	for (int simIndex = 0; simIndex < numSimulations; simIndex++)
	{
		sum += discountedPayoff(simIndex, numExerciseDates, basketValues, exerciseTimes);
	}

	double price = (sum / numSimulations);

	return price;
}

// Bermudan Helper Functions
void BermudanBasket::backwardInduction(int dateIndex, int numSimulations, Eigen::MatrixXd& basketValues, Eigen::MatrixXd& exerciseTimes)
{
	int numExerciseDates = exerciseDates.size();
    double currentExerciseTime = exerciseDates[dateIndex];
    Eigen::MatrixXd regressionMatrix(numSimulations, L);
    Eigen::MatrixXd regressionTarget(numSimulations, 1);

	// Calculate regression coefficients
    for (int simIndex = 0; simIndex < numSimulations; simIndex++)
    {
        double nextExerciseTime = exerciseTimes(simIndex, dateIndex + 1);
        int lastIndex = numExerciseDates - 1;
        while (nextExerciseTime != exerciseDates[lastIndex])
            lastIndex--;

        regressionTarget(simIndex, 0) = std::exp(-rates[0] * (exerciseTimes(simIndex, dateIndex + 1) - currentExerciseTime)) * std::max(basketValues(simIndex, lastIndex) - strike, 0.);

        for (int powerIndex = 0; powerIndex < L; powerIndex++)
            regressionMatrix(simIndex, powerIndex) = pow(basketValues(simIndex, dateIndex), powerIndex);
    }

	// Solve for regression coefficients
    Eigen::VectorXd regressionCoefficients(regressionMatrix.colPivHouseholderQr().solve(regressionTarget));

	// Determine optimal exercise strategy
    for (int simIndex = 0; simIndex < numSimulations; simIndex++)
    {
        double actualPayoff = std::max(basketValues(simIndex, dateIndex) - strike, 0.);
        double expectedPayoff = regressionMatrix.row(simIndex) * regressionCoefficients;

        exerciseTimes(simIndex, dateIndex) = (actualPayoff >= expectedPayoff) ? currentExerciseTime : exerciseTimes(simIndex, dateIndex + 1);
    }
}

double BermudanBasket::discountedPayoff(int simIndex, int numExerciseDates, Eigen::MatrixXd& basketValues, Eigen::MatrixXd& exerciseTimes)
{
	// Find exercise time at initial date
	double initialExerciseTime = exerciseTimes(simIndex, 0);
	int lastIndex = numExerciseDates - 1;

	while (initialExerciseTime != exerciseDates[lastIndex])
		lastIndex--;

	double discountedPayoff = std::exp(-rates[0] * initialExerciseTime) * std::max(basketValues(simIndex, lastIndex) - strike, 0.);

	paths_prices[simIndex] = discountedPayoff;
	return discountedPayoff;
}