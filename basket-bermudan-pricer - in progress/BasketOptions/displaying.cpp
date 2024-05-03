#include <iostream>
#include "displaying.h"

Display::Display(Priceable* pric) {}

void Display::display(Priceable* pric, bool quasiMC)
{
    std::pair<double, double> pricesRaw = pric->getPricesRaw();
    std::pair<double, double> pricesAntithetic;
    std::pair<double, double> pricesControlVar;

    if (!quasiMC) {
        pricesAntithetic = pric->getPricesAntithetic();
        pricesControlVar = pric->getPricesControlVar();
    }
    // no need for an else block for the quasiMC case

    std::string optionFeaturesName = pric->getOptionFeaturesName();

    // Print out the option names and prices
    std::cout << "Price of a European " << optionFeaturesName << " : " << pricesRaw.first << std::endl;
    std::cout << "Price of a Bermudan " << optionFeaturesName << " : " << pricesRaw.second << std::endl;

    if (!quasiMC) {
        // pricesAntithetic
        std::cout << "Price of a European " << optionFeaturesName << " with antithetic variable : " << pricesAntithetic.first << std::endl;
        std::cout << "Price of a Bermudan " << optionFeaturesName << " with antithetic variable : " << pricesAntithetic.second << std::endl;

        // pricesControlVar
        std::cout << "Price of a European " << optionFeaturesName << " with control variate : " << pricesControlVar.first << std::endl;
        std::cout << "Price of a Bermudan " << optionFeaturesName << " with control variate : " << pricesControlVar.second << std::endl;
    }

    std::cout << "\n" << std::endl;
}