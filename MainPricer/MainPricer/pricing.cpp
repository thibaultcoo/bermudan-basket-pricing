#include "pricing.h"

Priceable::Priceable(Instance inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption, std::string optionFeaturesName, bool quasiMC)
    : optionFeaturesName(optionFeaturesName) {
    computeMetrics(inst, EuropeanOption, BermudanOption, quasiMC);
}

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

std::pair<double, double> Priceable::computePricesRaw(Instance& inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption) {
    return { EuropeanOption->price(inst.getNbPaths()), BermudanOption->price(inst.getNbPaths()) };
}

std::pair<double, double> Priceable::computePricesAntithetic(Instance& inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption) {
    return { EuropeanOption->priceAntithetic(inst.getNbPaths()), BermudanOption->priceAntithetic(inst.getNbPaths()) };
}

std::pair<double, double> Priceable::computePricesControlVar(Instance& inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption) {
    return { EuropeanOption->priceControlVariate(inst.getNbPaths()), BermudanOption->priceControlVariate(inst.getNbPaths()) };
}

std::pair<double, double> Priceable::computeVariances(EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption) {
    return { EuropeanOption->variance(), BermudanOption->variance() };
}

std::pair<std::vector<double>, std::vector<double>> Priceable::computeConfInterval(EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption) {
    return { EuropeanOption->confidenceInterval(0.99), BermudanOption->confidenceInterval(0.99) };
}