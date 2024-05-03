#include "pricing.h"

Priceable::Priceable(Instance inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption, std::string optionFeaturesName, bool quasiMC)
    : optionFeaturesName(optionFeaturesName) {
    if (!quasiMC) {
        pricesRaw = computePricesRaw(inst, EuropeanOption, BermudanOption);
        pricesAntithetic = computePricesAntithetic(inst, EuropeanOption, BermudanOption);
        pricesControlVar = computePricesControlVar(inst, EuropeanOption, BermudanOption);
    }
    else {
        pricesRaw = computePricesRaw(inst, EuropeanOption, BermudanOption);
    }
}

std::pair<double, double> Priceable::getPricesRaw() const {
    return pricesRaw;
}

std::pair<double, double> Priceable::getPricesAntithetic() const {
    return pricesAntithetic;
}

std::pair<double, double> Priceable::getPricesControlVar() const {
    return pricesControlVar;
}

std::pair<double, double> Priceable::getPricesQuasiMC() const {
    return pricesQuasiMC;
}

std::string Priceable::getOptionFeaturesName() const {
    return optionFeaturesName;
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

std::pair<double, double> Priceable::computePricesQuasiMC(Instance& inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption) {
    return { EuropeanOption->priceQuasiMC(inst.getNbPaths()), BermudanOption->priceQuasiMC(inst.getNbPaths()) };
}