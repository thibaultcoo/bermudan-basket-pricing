#include "EuropeanBasket.h"
#include "BermudanBasket.h"
#include "instantiating.h"
#include "pricing.h"
#include "displaying.h"
#include <eigen-3.4.0/Eigen/Dense>

/*
Pricing of bermudan basket options - M203 Numerical finance project
Marchessaux François, Louis Faverjon, Collin Thibault, Piancastelli Robin
*/

int main()
{
    // Creating an instance of the complete algorithm
    Instance Inst;

    // Initializing the base parameter for the algorithm and the derivatives
    Inst.fillWithParameters();

    // Initializing the random sequence generators that will be used to diffuse the underlyings
    Inst.fillWithGenerators();

    // Diffusing the schemes and storing them in this full instance of an algorithm
    Inst.fillWithSchemes();

    // Initializing the European option objects
    EuropeanBasket* EuropeanEuler = new EuropeanBasket(Inst.getEuler(), Inst);
    EuropeanBasket* EuropeanMilstein = new EuropeanBasket(Inst.getMilstein(), Inst);
    EuropeanBasket* EuropeanEulerQuasiMC = new EuropeanBasket(Inst.getEulerQuasiMC(), Inst);
    EuropeanBasket* EuropeanMilsteinQuasiMC = new EuropeanBasket(Inst.getMilsteinQuasiMC(), Inst);

    // Initializing the Bermudan option objects
    BermudanBasket* BermudanEuler = new BermudanBasket(Inst.getEuler(), Inst);
    BermudanBasket* BermudanMilstein = new BermudanBasket(Inst.getMilstein(), Inst);
    BermudanBasket* BermudanEulerQuasiMC = new BermudanBasket(Inst.getEulerQuasiMC(), Inst);
    BermudanBasket* BermudanMilsteinQuasiMC = new BermudanBasket(Inst.getMilsteinQuasiMC(), Inst);

    // Finalizing the initialization the priceable objects
    std::cout << "Pricing European and Bermudan Basket options with Euler scheme ..." << std::endl;
    Priceable* Euler = new Priceable(Inst, EuropeanEuler, BermudanEuler, "Basket Call Euler", false);
    std::cout << "Pricing European and Bermudan Basket options with Milstein scheme ..." << std::endl;
    Priceable* Milstein = new Priceable(Inst, EuropeanMilstein, BermudanMilstein, "Basket Call Milstein", false);
    std::cout << "Pricing European and Bermudan Basket options with Euler scheme and Quasi Monte Carlo ..." << std::endl;
    Priceable* EulerQuasiMC = new Priceable(Inst, EuropeanEulerQuasiMC, BermudanEulerQuasiMC, "Basket Call Euler with Quasi Monte-Carlo", true);
    std::cout << "Pricing European and Bermudan Basket options with Milstein scheme and Quasi Monte Carlo ... \n\n\n" << std::endl;
    Priceable* MilsteinQuasiMC = new Priceable(Inst, EuropeanMilsteinQuasiMC, BermudanMilsteinQuasiMC, "Basket Call Milstein with Quasi Monte-Carlo", true);

    // Initializing the clean displaying prompter
    Display* EulerRes = new Display(Euler);
    Display* MilsteinRes = new Display(Milstein);
    Display* EulerQuasiMCRes = new Display(EulerQuasiMC);
    Display* MilsteinQuasiMCRes = new Display(MilsteinQuasiMC);

    // Displaying the pricing logs and results
    EulerRes->display(Euler, false);
    MilsteinRes->display(Milstein, false);
    EulerQuasiMCRes->display(EulerQuasiMC, true);
    MilsteinQuasiMCRes->display(MilsteinQuasiMC, true);

    // Displaying our instance initial parameters
    paramsDisplay* params = new paramsDisplay(Inst);
    params->display(Inst);

    return 0;
}
