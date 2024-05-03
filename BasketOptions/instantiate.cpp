#include <iostream>
#include "LinearCongruential.h"
#include "VanDerCorput.h"
#include "EcuyerCombined.h"
#include "NormalBoxMuller.h"
#include "NormalCLT.h"
#include "NormalInvCdf.h"
#include "BSEulernD.h"
#include "BSMilsteinnD.h"
#include "instantiate.h"
#include <eigen-3.4.0/Eigen/Dense>

// Instance object constructor
instance::instance()
{
    // Constructor body left empty since constants are initialized later within several functions.
}

void instance::fillWithParameters()
{
    i_bermudans = { 1.0 / 12.0, 1.0 / 6.0, 3.0 / 12.0 };
    i_spots = { 30.0, 40.0, 5.0 };
    i_rates = { 0.05, 0.05, 0.05 };
    i_weights = { 0.25, 0.25, 0.5 };

    i_matu = 3.0 / 12.0;
    i_strike = 4.0;
    i_numAssets = i_spots.size();
    i_nbPaths = 1000;
    i_nbPathsQuasiMC = 200;

    // Redefine variances matrix
    i_variances.resize(3, 3);
    i_variances << 0.06, 0.04, 0.03,
        0.04, 0.07, 0.045,
        0.03, 0.045, 0.05;
}

void instance::fillWithGenerators()
{
    Unif = new EcuyerCombined();
    QuasiMC = new VanDerCorput();
    Normal_Unif = new NormalCLT(0, 1, Unif);
    NormalQuasiMC = new NormalInvCdf(0, 1, QuasiMC);
}

void instance::fillWithSchemes()
{
    Euler = new BSEuler(Normal_Unif, i_spots, i_rates, i_variances, i_numAssets);
    EulerQuasiMC = new BSEuler(NormalQuasiMC, i_spots, i_rates, i_variances, i_numAssets);
    Milstein = new BSMilstein(Normal_Unif, i_spots, i_rates, i_variances, i_numAssets);
    MilsteinQuasiMC = new BSMilstein(NormalQuasiMC, i_spots, i_rates, i_variances, i_numAssets);
}