#pragma once

#include <iostream>
#include "LinearCongruential.h"
#include "VanDerCorput.h"
#include "EcuyerCombined.h"
#include "Normal.h"
#include "BSEuler.h"
#include "BSMilstein.h"
#include <eigen-3.4.0/Eigen/Dense>

class Instance
{
public:
    Instance();
    void fillWithParameters();
    void fillWithGenerators();
    void fillWithSchemes();

    // Getter methods for member variables and objectss
    const std::vector<double>& getBermudans() const { return i_bermudans; }
    const std::vector<double>& getSpots() const { return i_spots; }
    const std::vector<double>& getRates() const { return i_rates; }
    const std::vector<double>& getWeights() const { return i_weights; }
    double getMatu() const { return i_matu; }
    double getStrike() const { return i_strike; }
    double getNumAssets() const { return i_numAssets; }
    int getNbPaths() const { return i_nbPaths; }
    int getNbPathsQuasiMC() const { return i_nbPathsQuasiMC; }
    const Eigen::MatrixXd& getVariances() const { return i_variances; }
    BSEuler* getEuler() const { return Euler; }
    BSEuler* getEulerQuasiMC() const { return EulerQuasiMC; }
    BSMilstein* getMilstein() const { return Milstein; }
    BSMilstein* getMilsteinQuasiMC() const { return MilsteinQuasiMC; }
    EcuyerCombined* getUnif() const { return Unif; }
    VanDerCorput* getQuasiMC() const { return QuasiMC; }
    NormalCLT* getNormalUnif() const { return Normal_Unif; }
    NormalInvCdf* getNormalQuasiMC() const { return NormalQuasiMC; }

private:
    std::vector<double> i_bermudans;
    std::vector<double> i_spots;
    std::vector<double> i_rates;
    std::vector<double> i_weights;

    double i_matu;
    double i_strike;
    double i_numAssets;
    int i_nbPaths;
    int i_nbPathsQuasiMC;

    Eigen::MatrixXd i_variances;

    // Pointers to schemes
    BSEuler* Euler;
    BSEuler* EulerQuasiMC;
    BSMilstein* Milstein;
    BSMilstein* MilsteinQuasiMC;

    // Pointers to generators
    EcuyerCombined* Unif;
    VanDerCorput* QuasiMC;
    NormalCLT* Normal_Unif;
    NormalInvCdf* NormalQuasiMC;
};
