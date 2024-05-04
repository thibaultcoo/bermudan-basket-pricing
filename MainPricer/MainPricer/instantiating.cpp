#include <iostream>
#include "LinearCongruential.h"
#include "VanDerCorput.h"
#include "EcuyerCombined.h"
#include "Normal.h"
#include "BSEuler.h"
#include "BSMilstein.h"
#include "instantiating.h"
#include <eigen-3.4.0/Eigen/Dense>

// Instance object constructor
Instance::Instance() {}

void Instance::fillWithParameters()
{
    i_bermudans = { 1.0 / 12.0, 1.0 / 6.0, 3.0 / 12.0 };
    i_spots = { 5.0, 10.0, 20.0 };
    i_rates = { 0.05, 0.05, 0.05 };
    i_weights = { 0.25, 0.25, 0.5 };

    i_matu = 3.0 / 12.0;
    i_strike = 10.0;
    i_numAssets = i_spots.size();
    i_nbPaths = 1000;
    i_nbPathsQuasiMC = 500;

    // Redefine variances matrix
    i_variances.resize(3, 3);
    i_variances << 0.06, 0.04, 0.03,
        0.04, 0.07, 0.045,
        0.03, 0.045, 0.05;
}

void Instance::fillWithGenerators()
{
    Unif = new EcuyerCombined();
    QuasiMC = new VanDerCorput();
    Normal_Unif = new NormalCLT(0, 1, Unif);
    NormalQuasiMC = new NormalInvCdf(0, 1, QuasiMC);
}

void Instance::fillWithSchemes()
{
    Euler = new BSEuler(Normal_Unif, i_spots, i_rates, i_variances, i_numAssets);
    EulerQuasiMC = new BSEuler(NormalQuasiMC, i_spots, i_rates, i_variances, i_numAssets);
    Milstein = new BSMilstein(Normal_Unif, i_spots, i_rates, i_variances, i_numAssets);
    MilsteinQuasiMC = new BSMilstein(NormalQuasiMC, i_spots, i_rates, i_variances, i_numAssets);
}