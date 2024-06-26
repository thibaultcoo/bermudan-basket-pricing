#include <iostream>
#include "LinearCongruential.h"
#include "VanDerCorput.h"
#include "EcuyerCombined.h"
#include "Normal.h"
#include "BSEuler.h"
#include "BSMilstein.h"
#include "instantiating.h"
#include <eigen-3.4.0/Eigen/Dense>

// Instance object constructor: it is meant to hold the main parameters, sequences and schemes of the resolution method
Instance::Instance() {}

// This function simply fills the instance with the base parameters for the option pricing: (Prof.) can freely modify any parameters here
void Instance::fillWithParameters()
{
    i_bermudans = { 1.0 / 12.0, 1.0 / 6.0, 3.0 / 12.0 };
    i_spots = { 5.0, 10.0, 20.0 };
    i_rates = { 0.01, 0.01, 0.01 };
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

// This function stores the random number generators in the instance object for further use (in the Euler/Milstein stochastic component of the path generation)
void Instance::fillWithGenerators()
{
    Unif = new EcuyerCombined();
    QuasiMC = new VanDerCorput();
    Normal_Unif = new NormalCLT(0, 1, Unif);
    NormalQuasiMC = new NormalInvCdf(0, 1, QuasiMC);
}

// This function builds the different schemes of trajectory diffusion, considering the previously generated random sequences
void Instance::fillWithSchemes()
{
    Euler = new BSEuler(Normal_Unif, i_spots, i_rates, i_variances, i_numAssets);
    EulerQuasiMC = new BSEuler(NormalQuasiMC, i_spots, i_rates, i_variances, i_numAssets);
    Milstein = new BSMilstein(Normal_Unif, i_spots, i_rates, i_variances, i_numAssets);
    MilsteinQuasiMC = new BSMilstein(NormalQuasiMC, i_spots, i_rates, i_variances, i_numAssets);
}
