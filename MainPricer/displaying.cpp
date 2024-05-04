#include <iostream>
#include <iomanip>
#include <vector>
#include "displaying.h"

// This element will display all the pricing results on a windows terminal for the user to see
Display::Display(Priceable* pric) {}

// This is the function that prompts out the results
void Display::display(Priceable* pric, bool quasiMC)
{
    // We initialize our variables first (storing the QuasiMC result in Raw for conveniency, it has no impact)
    std::pair<double, double> pricesRaw;
    std::pair<double, double> pricesAntithetic;
    std::pair<double, double> pricesControlVar;

    std::pair<double, double> variancesRaw;
    std::pair<double, double> variancesAntithetic;
    std::pair<double, double> variancesControlVar;

    std::pair<std::vector<double>, std::vector<double>> confIntervalsRaw;
    std::pair<std::vector<double>, std::vector<double>> confIntervalsAntithetic;
    std::pair<std::vector<double>, std::vector<double>> confIntervalsControlVar;

    // In that case, we also consider Monte-Carlo variations, so we include their pricing results
    if (!quasiMC) {
        pricesRaw = pric->getPricesRaw();
        pricesAntithetic = pric->getPricesAntithetic();
        pricesControlVar = pric->getPricesControlVar();

        variancesRaw = pric->getVariancesRaw();
        variancesAntithetic = pric->getVariancesAntithetic();
        variancesControlVar = pric->getVariancesControlVar();

        confIntervalsRaw = pric->getConfIntervalRaw();
        confIntervalsAntithetic = pric->getConfIntervalAntithetic();
        confIntervalsControlVar = pric->getConfIntervalControlVar();
    }
    else {
        pricesRaw = pric->getPricesQuasiMC();
        variancesRaw = pric->getVariancesQuasiMC();
        confIntervalsRaw = pric->getConfIntervalQuasiMC();
    }

    // The option name is also retrieved to increase readability
    std::string optionFeaturesName = pric->getOptionFeaturesName();

    // Setup formatting for the output.
    const int nameWidth = 75;
    const int priceWidth = 12;
    const int varianceWidth = 12;
    const int ciWidth = 20;
    const char separator = ' ';

    // Priting the header for the table
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << "Option Type";
    std::cout << std::right << std::setw(priceWidth) << std::setfill(separator) << "Price";
    std::cout << std::right << std::setw(varianceWidth) << std::setfill(separator) << "Variance";
    std::cout << std::right << std::setw(ciWidth) << std::setfill(separator) << "IC" << std::endl;
    std::cout << std::setfill('-') << std::setw(nameWidth + priceWidth + varianceWidth + ciWidth) << "-" << std::endl;

    // Displaying details for each option
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

// This is a helper function to help structure and format the table
void Display::displayDetails(const std::string& description, double price, double variance, const std::vector<double>& confInterval, int nameWidth, int priceWidth, int varianceWidth, int ciWidth)
{
    // Printing the description
    std::cout << std::left << std::setw(nameWidth) << std::setfill(' ') << description;

    // Printing the price with proper formatting
    std::cout << std::right << std::setw(priceWidth) << std::setfill(' ') << price;

    // Printing the variance with proper formatting
    std::cout << std::right << std::setw(varianceWidth) << std::setfill(' ') << variance;

    // Preparing and printing the confidence interval string with proper formatting
    std::ostringstream ciStream;
    ciStream << std::fixed << std::setprecision(2) << "[" << confInterval[0] << ", " << confInterval[1] << "]";
    std::string ciStr = ciStream.str();
    std::cout << std::right << std::setw(ciWidth) << std::setfill(' ') << ciStr;

    std::cout << std::endl;
}

// This is another function called to prompt out the initial parameters used for the pricing: encompasses the derivatives market parameters as well as the MC simulation hyperparameters
void paramsDisplay::display(const Instance& inst) {
    std::cout << "\n" << std::endl;
    std::cout << "Model Resolution Parameters:" << std::endl;
    std::cout << "--------------------------------------" << std::endl;

    // Displaying the parameters
    std::cout << "Maturity (months): " << inst.getMatu() * 12 << std::endl;
    std::cout << "Strike price: " << inst.getStrike() << std::endl;
    std::cout << "Number of Assets: " << inst.getNumAssets() << std::endl;
    std::cout << "Number of Paths: " << inst.getNbPaths() << std::endl;
    std::cout << "Number of Paths for Quasi-MC: " << inst.getNbPathsQuasiMC() << std::endl;
    std::cout << "Rates: ";
    for (auto rate : inst.getRates()) std::cout << rate << " ";
    std::cout << std::endl;

    std::cout << "Spots: ";
    for (auto spot : inst.getSpots()) std::cout << spot << " ";
    std::cout << std::endl;

    std::cout << "Weights: ";
    for (auto weight : inst.getWeights()) std::cout << weight << " ";
    std::cout << std::endl;

    std::cout << "Bermudan Exercise Times (months): ";
    for (auto time : inst.getBermudans()) std::cout << time * 12 << " ";
    std::cout << std::endl;

    std::cout << "Variances Matrix:" << std::endl;
    for (int i = 0; i < inst.getVariances().rows(); ++i) {
        for (int j = 0; j < inst.getVariances().cols(); ++j) {
            std::cout << std::setw(5) << std::fixed << std::setprecision(3) << inst.getVariances()(i, j) << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "--------------------------------------" << std::endl;
}