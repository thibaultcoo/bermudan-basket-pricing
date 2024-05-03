#include <iostream>
#include "EuropeanBasketCall.h"
#include "BermudanBasketCall.h"
#include "instantiate.h"
#include <eigen-3.4.0/Eigen/Dense>

int main() {

    // Creating an instance of the complete algorithm.
    instance inst;

    // Initializing the base parameter for the algorithm and the derivatives.
    inst.fillWithParameters();

    // Initializing the random sequence generators that will be used to diffuse the underlyings.
    inst.fillWithGenerators();

    // Diffusing the schemes and storing them in this full instance of an algorithm.
    inst.fillWithSchemes();
   
    // Initializing the options.
    European* EuropeanEuler = new European(inst.getEuler(), inst.getStrike(), inst.getRates(), inst.getMatu(), inst.getWeights(), inst.getSpots(), inst.getVariances());
    European* EuropeanEulerQuasiMC = new European(inst.getEulerQuasiMC(), inst.getStrike(), inst.getRates(), inst.getMatu(), inst.getWeights(), inst.getSpots(), inst.getVariances());
    European* EuropeanMilstein = new European(inst.getMilstein(), inst.getStrike(), inst.getRates(), inst.getMatu(), inst.getWeights(), inst.getSpots(), inst.getVariances());
    European* EuropeanMilsteinQuasiMC = new European(inst.getMilsteinQuasiMC(), inst.getStrike(), inst.getRates(), inst.getMatu(), inst.getWeights(), inst.getSpots(), inst.getVariances());

    Bermudan* BermudanEuler = new Bermudan(inst.getEuler(), inst.getStrike(), inst.getRates(), inst.getMatu(), inst.getWeights(), inst.getBermudans(), inst.getSpots(), inst.getVariances());
    Bermudan* BermudanEulerQuasiMC = new Bermudan(inst.getEulerQuasiMC(), inst.getStrike(), inst.getRates(), inst.getMatu(), inst.getWeights(), inst.getBermudans(), inst.getSpots(), inst.getVariances());
    Bermudan* BermudanMilstein = new Bermudan(inst.getMilstein(), inst.getStrike(), inst.getRates(), inst.getMatu(), inst.getWeights(), inst.getBermudans(), inst.getSpots(), inst.getVariances());
    Bermudan* BermudanMilsteinQuasiMC = new Bermudan(inst.getMilsteinQuasiMC(), inst.getStrike(), inst.getRates(), inst.getMatu(), inst.getWeights(), inst.getBermudans(), inst.getSpots(), inst.getVariances());

    // below is the pricing and direct output of the values (we will split the pricing/displaying and aggregate the displaying nicely by introducing result objects)
    std::cout << "Price European Basket Call Euler simple : " << EuropeanEuler->ComputePrice(inst.getNbPaths()) << std::endl;
    std::cout << " - Variance : " << EuropeanEuler->calculate_variance() << std::endl;
    std::cout << " - IC : [" << EuropeanEuler->calculate_ConfidenceInterval()[0] << ", ";
    std::cout << EuropeanEuler->calculate_ConfidenceInterval()[1] << "]" << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price European Basket Call Milstein simple : " << EuropeanMilstein->ComputePrice(inst.getNbPaths()) << std::endl;
    std::cout << " - Variance : " << EuropeanMilstein->calculate_variance() << std::endl;
    std::cout << " - IC : [" << EuropeanMilstein->calculate_ConfidenceInterval()[0] << ", ";
    std::cout << EuropeanMilstein->calculate_ConfidenceInterval()[1] << "]" << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price Bermudan Basket Call Euler simple : " << BermudanEuler->ComputePrice(inst.getNbPaths()) << std::endl;
    std::cout << " - Variance : " << BermudanEuler->calculate_variance() << std::endl;
    std::cout << " - IC : [" << BermudanEuler->calculate_ConfidenceInterval()[0] << ", ";
    std::cout << BermudanEuler->calculate_ConfidenceInterval()[1] << "]" << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price Bermudan Basket Call Milstein simple : " << BermudanMilstein->ComputePrice(inst.getNbPaths()) << std::endl;
    std::cout << " - Variance : " << BermudanMilstein->calculate_variance() << std::endl;
    std::cout << " - IC : [" << BermudanMilstein->calculate_ConfidenceInterval()[0] << ", ";
    std::cout << BermudanMilstein->calculate_ConfidenceInterval()[1] << "]" << std::endl;
    std::cout << "\n" << std::endl;

    return 0;
}