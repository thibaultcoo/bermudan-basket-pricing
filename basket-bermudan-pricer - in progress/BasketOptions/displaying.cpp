#include <iostream>
#include <iomanip>
#include "displaying.h"

Display::Display(Priceable* pric) : priceable(pric) {}

void Display::display(Priceable* pric, bool quasiMC)
{
    std::pair<double, double> pricesRaw = pric->getPricesRaw();
    std::pair<double, double> pricesAntithetic;
    std::pair<double, double> pricesControlVar;

    std::pair<double, double> variancesRaw = pric->getVariancesRaw();
    std::pair<double, double> variancesAntithetic;
    std::pair<double, double> variancesControlVar;

    std::pair<std::vector<double>, std::vector<double>> confIntervalsRaw = pric->getConfIntervalRaw();
    std::pair<std::vector<double>, std::vector<double>> confIntervalsAntithetic;
    std::pair<std::vector<double>, std::vector<double>> confIntervalsControlVar;

    if (!quasiMC) {
        pricesAntithetic = pric->getPricesAntithetic();
        pricesControlVar = pric->getPricesControlVar();

        variancesAntithetic = pric->getVariancesAntithetic();
        variancesControlVar = pric->getVariancesControlVar();

        confIntervalsAntithetic = pric->getConfIntervalAntithetic();
        confIntervalsControlVar = pric->getConfIntervalControlVar();
    }

    std::string optionFeaturesName = pric->getOptionFeaturesName();

    // Setup formatting for the output
    const int nameWidth = 75;
    const int priceWidth = 12;
    const int varianceWidth = 12;
    const int ciWidth = 20;
    const char separator = ' ';

    // Print a header
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << "Option Type";
    std::cout << std::right << std::setw(priceWidth) << std::setfill(separator) << "Price";
    std::cout << std::right << std::setw(varianceWidth) << std::setfill(separator) << "Variance";
    std::cout << std::right << std::setw(ciWidth) << std::setfill(separator) << "IC" << std::endl;

    std::cout << std::setfill('-') << std::setw(nameWidth + priceWidth + varianceWidth + ciWidth) << "-" << std::endl;

    // Display details for each option
    displayDetails("European " + optionFeaturesName + ":", pricesRaw.first, variancesRaw.first, confIntervalsRaw.first, nameWidth, priceWidth, varianceWidth, ciWidth);
    displayDetails("Bermudan " + optionFeaturesName + ":", pricesRaw.second, variancesRaw.second, confIntervalsRaw.second, nameWidth, priceWidth, varianceWidth, ciWidth);

    if (!quasiMC) {
        displayDetails("European " + optionFeaturesName + " with antithetic variable:", pricesAntithetic.first, variancesAntithetic.first, confIntervalsAntithetic.first, nameWidth, priceWidth, varianceWidth, ciWidth);
        displayDetails("Bermudan " + optionFeaturesName + " with antithetic variable:", pricesAntithetic.second, variancesAntithetic.second, confIntervalsAntithetic.second, nameWidth, priceWidth, varianceWidth, ciWidth);

        displayDetails("European " + optionFeaturesName + " with control variate:", pricesControlVar.first, variancesControlVar.first, confIntervalsControlVar.first, nameWidth, priceWidth, varianceWidth, ciWidth);
        displayDetails("Bermudan " + optionFeaturesName + " with control variate:", pricesControlVar.second, variancesControlVar.second, confIntervalsControlVar.second, nameWidth, priceWidth, varianceWidth, ciWidth);
    }

    std::cout << std::setfill('-') << std::setw(nameWidth + priceWidth + varianceWidth + ciWidth) << "-" << std::endl;
    std::cout << "\n";
}

void Display::displayDetails(const std::string& description, double price, double variance, const std::vector<double>& confInterval, int nameWidth, int priceWidth, int varianceWidth, int ciWidth) {
    // Print the description
    std::cout << std::left << std::setw(nameWidth) << std::setfill(' ') << description;

    // Print the price with proper formatting
    std::cout << std::right << std::setw(priceWidth) << std::setfill(' ') << price;

    // Print the variance with proper formatting
    std::cout << std::right << std::setw(varianceWidth) << std::setfill(' ') << variance;

    // Prepare and print the confidence interval string with proper formatting
    std::ostringstream ciStream;
    ciStream << std::fixed << std::setprecision(2) << "[" << confInterval[0] << ", " << confInterval[1] << "]";
    std::string ciStr = ciStream.str();
    std::cout << std::right << std::setw(ciWidth) << std::setfill(' ') << ciStr;

    std::cout << std::endl;
}