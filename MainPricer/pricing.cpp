#include "pricing.h"

// Initialization of the priceable object that is meant to build and price some derivatives, alongside their MC variance and confidence interval
Priceable::Priceable(Instance inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption, std::string optionFeaturesName, bool quasiMC)
    : optionFeaturesName(optionFeaturesName) {
    computeMetrics(inst, EuropeanOption, BermudanOption, quasiMC);
}

// Function used to encompass the computation of variances and confidence intervals for all the required pricing variations for a given structure
void Priceable::computeMetrics(Instance& inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption, bool quasiMC) {
    pricesRaw = computePricesRaw(inst, EuropeanOption, BermudanOption);
    variancesRaw = computeVariances(EuropeanOption, BermudanOption);
    confIntervalRaw = computeConfInterval(EuropeanOption, BermudanOption);

    if (!quasiMC) {
        pricesAntithetic = computePricesAntithetic(inst, EuropeanOption, BermudanOption);
        variancesAntithetic = computeVariances(EuropeanOption, BermudanOption);
        confIntervalAntithetic = computeConfInterval(EuropeanOption, BermudanOption);

        pricesControlVar = computePricesControlVar(inst, EuropeanOption, BermudanOption);
        variancesControlVar = computeVariances(EuropeanOption, BermudanOption);
        confIntervalControlVar = computeConfInterval(EuropeanOption, BermudanOption);
    }
}

// Function that returns the raw prices of both the European and Bermudan options as a pair
std::pair<double, double> Priceable::computePricesRaw(Instance& inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption) {
    return { EuropeanOption->price(inst.getNbPaths()), BermudanOption->price(inst.getNbPaths()) };
}

// Function that returns the prices of both the European and Bermudan options as a pair considering adding an antithetic variable for variance reduction purposes
std::pair<double, double> Priceable::computePricesAntithetic(Instance& inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption) {
    return { EuropeanOption->priceAntithetic(inst.getNbPaths()), BermudanOption->priceAntithetic(inst.getNbPaths()) };
}

// Function that returns the prices of both the European and Bermudan options as a pair considering adding a control variate for variance reduction purposes
std::pair<double, double> Priceable::computePricesControlVar(Instance& inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption) {
    return { EuropeanOption->priceControlVariate(inst.getNbPaths()), BermudanOption->priceControlVariate(inst.getNbPaths()) };
}

// Function that returns the pair of variances for both options
std::pair<double, double> Priceable::computeVariances(EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption) {
    return { EuropeanOption->variance(), BermudanOption->variance() };
}

// Function that returns the pair of confidence intervals for both options
std::pair<std::vector<double>, std::vector<double>> Priceable::computeConfInterval(EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption) {
    return { EuropeanOption->confidenceInterval(0.99), BermudanOption->confidenceInterval(0.99) };
}