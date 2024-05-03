#include <iostream>
#include "EUBasketCall.h"
#include "BermudanBasket.h"
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


    EUBasketCall* EuropeanEuler = new EUBasketCall(inst.getEuler(), inst.getStrike(), inst.getRates(), inst.getMatu(), inst.getWeights(), inst.getSpots(), inst.getVariances());
    EUBasketCall* EuropeanEulerQuasiMC = new EUBasketCall(inst.getEulerQuasiMC(), inst.getStrike(), inst.getRates(), inst.getMatu(), inst.getWeights(), inst.getSpots(), inst.getVariances());
    EUBasketCall* EuropeanMilstein = new EUBasketCall(inst.getMilstein(), inst.getStrike(), inst.getRates(), inst.getMatu(), inst.getWeights(), inst.getSpots(), inst.getVariances());
    EUBasketCall* EuropeanMilsteinQuasiMC = new EUBasketCall(inst.getMilsteinQuasiMC(), inst.getStrike(), inst.getRates(), inst.getMatu(), inst.getWeights(), inst.getSpots(), inst.getVariances());

    BermudanBasket* BermudanEuler = new BermudanBasket(inst.getEuler(), inst.getStrike(), inst.getRates(), inst.getMatu(), inst.getWeights(), inst.getBermudans(), inst.getSpots(), inst.getVariances());
    BermudanBasket* BermudanEulerQuasiMC = new BermudanBasket(inst.getEulerQuasiMC(), inst.getStrike(), inst.getRates(), inst.getMatu(), inst.getWeights(), inst.getBermudans(), inst.getSpots(), inst.getVariances());
    BermudanBasket* BermudanMilstein = new BermudanBasket(inst.getMilstein(), inst.getStrike(), inst.getRates(), inst.getMatu(), inst.getWeights(), inst.getBermudans(), inst.getSpots(), inst.getVariances());
    BermudanBasket* BermudanMilsteinQuasiMC = new BermudanBasket(inst.getMilsteinQuasiMC(), inst.getStrike(), inst.getRates(), inst.getMatu(), inst.getWeights(), inst.getBermudans(), inst.getSpots(), inst.getVariances());

    // below is the pricing and direct output of the values (we will split the pricing/displaying and aggregate the displaying nicely by introducing result objects)
    std::cout << "Price European Basket Call Euler simple : " << EuropeanEuler->priceAntithetic(inst.getNbPaths()) << std::endl;
    std::cout << " - Variance : " << EuropeanEuler->variance() << std::endl;
    std::cout << " - IC : [" << EuropeanEuler->confidenceInterval()[0] << ", ";
    std::cout << EuropeanEuler->confidenceInterval()[1] << "]" << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price European Basket Call Milstein simple : " << EuropeanMilstein->price(inst.getNbPaths()) << std::endl;
    std::cout << " - Variance : " << EuropeanMilstein->variance() << std::endl;
    std::cout << " - IC : [" << EuropeanMilstein->confidenceInterval()[0] << ", ";
    std::cout << EuropeanMilstein->confidenceInterval()[1] << "]" << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price Bermudan Basket Call Euler simple : " << BermudanEuler->priceControlVariate(inst.getNbPaths()) << std::endl;
    std::cout << " - Variance : " << BermudanEuler->variance() << std::endl;
    std::cout << " - IC : [" << BermudanEuler->confidenceInterval()[0] << ", ";
    std::cout << BermudanEuler->confidenceInterval()[1] << "]" << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << "Price Bermudan Basket Call Milstein simple : " << BermudanMilstein->price(inst.getNbPaths()) << std::endl;
    std::cout << " - Variance : " << BermudanMilstein->variance() << std::endl;
    std::cout << " - IC : [" << BermudanMilstein->confidenceInterval()[0] << ", ";
    std::cout << BermudanMilstein->confidenceInterval()[1] << "]" << std::endl;
    std::cout << "\n" << std::endl;

    return 0;
}