#pragma once

#include <iostream>
#include <utility>
#include <vector>
#include "instantiating.h"
#include "EuropeanBasket.h"
#include "BermudanBasket.h"

class Priceable {
public:
    // Constructor that initializes the Priceable object.
    Priceable(Instance inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption, std::string optionFeaturesName, bool quasiMC);

    // Getter methods for prices
    std::pair<double, double> getPricesRaw() const { return pricesRaw; }
    std::pair<double, double> getPricesAntithetic() const { return pricesAntithetic; }
    std::pair<double, double> getPricesControlVar() const { return pricesControlVar; }

    // Getter methods for variance and confidence intervals for each price type
    std::pair<double, double> getVariancesRaw() const { return variancesRaw; }
    std::pair<double, double> getVariancesAntithetic() const { return variancesAntithetic; }
    std::pair<double, double> getVariancesControlVar() const { return variancesControlVar; }

    std::pair<std::vector<double>, std::vector<double>> getConfIntervalRaw() const { return confIntervalRaw; }
    std::pair<std::vector<double>, std::vector<double>> getConfIntervalAntithetic() const { return confIntervalAntithetic; }
    std::pair<std::vector<double>, std::vector<double>> getConfIntervalControlVar() const { return confIntervalControlVar; }

    // Getter for option features name
    std::string getOptionFeaturesName() const { return optionFeaturesName; }

private:
    // Member variables to store prices
    std::pair<double, double> pricesRaw, pricesAntithetic, pricesControlVar;
    std::string optionFeaturesName;

    // Variables to store variances and confidence intervals for each price type
    std::pair<double, double> variancesRaw, variancesAntithetic, variancesControlVar;
    std::pair<std::vector<double>, std::vector<double>> confIntervalRaw, confIntervalAntithetic, confIntervalControlVar;

    // Compute functions for prices
    std::pair<double, double> computePricesRaw(Instance& inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption);
    std::pair<double, double> computePricesAntithetic(Instance& inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption);
    std::pair<double, double> computePricesControlVar(Instance& inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption);

    // Methods to compute variances and confidence intervals
    std::pair<double, double> computeVariances(EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption);
    std::pair<std::vector<double>, std::vector<double>> computeConfInterval(EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption);

    // A method to compute all metrics including variances and confidence intervals
    void computeMetrics(Instance& inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption, bool quasiMC);
};