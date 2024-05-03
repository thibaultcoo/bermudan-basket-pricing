#pragma once

#include <iostream>
#include <utility>
#include "instantiating.h"
#include "EuropeanBasket.h"
#include "BermudanBasket.h"

class Priceable
{
public:
    Priceable(Instance inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption, std::string optionFeaturesName, bool quasiMC);

    // Getter methods to retrieve the prices and option names
    std::pair<double, double> getPricesRaw() const;
    std::pair<double, double> getPricesAntithetic() const;
    std::pair<double, double> getPricesControlVar() const;
    std::pair<double, double> getPricesQuasiMC() const;
    std::string getOptionFeaturesName() const;

private:
    std::pair<double, double> pricesRaw;
    std::pair<double, double> pricesAntithetic;
    std::pair<double, double> pricesControlVar;
    std::pair<double, double> pricesQuasiMC;
    std::string optionFeaturesName;

    std::pair<double, double> computePricesRaw(Instance& inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption);
    std::pair<double, double> computePricesAntithetic(Instance& inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption);
    std::pair<double, double> computePricesControlVar(Instance& inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption);
    std::pair<double, double> computePricesQuasiMC(Instance& inst, EuropeanBasket* EuropeanOption, BermudanBasket* BermudanOption);
};