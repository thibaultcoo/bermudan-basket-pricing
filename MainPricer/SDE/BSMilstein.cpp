#include "pch.h"
#include "BSMilstein.h"
#include <cmath>
#include <algorithm>


BSMilstein::BSMilstein() { }

BSMilstein::BSMilstein(RandomGenerator* generator, std::vector<double> spots, std::vector<double> rates, Eigen::MatrixXd corrMatrix, int dimension) : BlackScholes(generator, spots, rates, corrMatrix, dimension)
{
    // Check if the matrix is singular
    if (corrMatrix.determinant() == 0)
    {
        Eigen::EigenSolver<Eigen::MatrixXd> es(corrMatrix);
        choleskyDecomposition = es.pseudoEigenvectors() * es.pseudoEigenvalueMatrix();
    }
    else
    {
        // Cholesky decomposition for non-singular matrices
        choleskyDecomposition = corrMatrix.llt().matrixL();
    }
}

void BSMilstein::Simulate(double start_time, double end_time, size_t nb_steps)
{
    double timeStep = (end_time - start_time) / nb_steps;
    Eigen::VectorXd currentValues = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(spots.data(), spots.size());

    // Initialize Trajectories
    for (int i = 0; i < dim; i++)
    {
        trajectory[i] = new SingleTrajectory(start_time, end_time, nb_steps);
        trajectory[i]->AddValue(currentValues[i]);
    }

    for (int i = 0; i < nb_steps; i++)
    {
        // BM Component
        std::vector<double> dW = gen->GenerateVector(dim);
        Eigen::VectorXd mapped_dW = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(dW.data(), dW.size())* sqrt(timeStep);

        // Drift Component
        Eigen::VectorXd mapped_drift = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(rates.data(), rates.size());
        Eigen::MatrixXd choleskySquare = choleskyDecomposition * choleskyDecomposition.transpose();
        Eigen::VectorXd milsteinCorrection = 0.5 * choleskySquare.diagonal();

        // Diffusion
        Eigen::VectorXd Z = choleskyDecomposition * mapped_dW;
        Eigen::VectorXd X = (mapped_drift - milsteinCorrection) * timeStep + Z + 0.5 * (Z.cwiseProduct(Z));
        currentValues += currentValues.cwiseProduct(X);
        for (int j = 0; j < dim; j++)
            trajectory[j]->AddValue(currentValues[j]);
    }
}

void BSMilstein::SimulateAntithetic(double start_time, double end_time, size_t nb_steps)
{
    double timeStep = (end_time - start_time) / nb_steps;
    Eigen::VectorXd currentValues = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(spots.data(), spots.size());
    Eigen::VectorXd currentValuesAnti = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(spots.data(), spots.size());

    // Initialize Trajectories
    for (int i = 0; i < dim; i++)
    {
        trajectory[i] = new SingleTrajectory(start_time, end_time, nb_steps);
        trajectory[i]->AddValue(currentValues[i]);

        trajectoryAntithetic[i] = new SingleTrajectory(start_time, end_time, nb_steps);
        trajectoryAntithetic[i]->AddValue(currentValuesAnti[i]);
    }

    for (int i = 0; i < nb_steps; i++)
    {
        // BM Component
        std::vector<double> dW = gen->GenerateVector(dim);
        Eigen::VectorXd mapped_dW = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(dW.data(), dW.size())* sqrt(timeStep);
        Eigen::VectorXd mapped_dW_anti = -mapped_dW;

        // Drift Component
        Eigen::VectorXd mapped_drift = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(rates.data(), rates.size());
        Eigen::MatrixXd choleskySquare = choleskyDecomposition * choleskyDecomposition.transpose();
        Eigen::VectorXd milsteinCorrection = 0.5 * choleskySquare.diagonal();

        // Diffusion
        Eigen::VectorXd Z = choleskyDecomposition * mapped_dW;
        Eigen::VectorXd Z_anti = choleskyDecomposition * mapped_dW_anti;
        Eigen::VectorXd X = (mapped_drift - milsteinCorrection) * timeStep + Z + 0.5 * (Z.cwiseProduct(Z));
        Eigen::VectorXd X_anti = (mapped_drift - milsteinCorrection) * timeStep + Z_anti + 0.5 * (Z_anti.cwiseProduct(Z_anti));
        currentValues += currentValues.cwiseProduct(X);
        currentValuesAnti += currentValuesAnti.cwiseProduct(X);
        for (int j = 0; j < dim; j++)
        {
            trajectory[j]->AddValue(currentValues[j]);
            trajectoryAntithetic[j]->AddValue(currentValuesAnti[j]);
        }
    }
}

void BSMilstein::SimulateQuasiMC(double start_time, double end_time, size_t nb_steps, myLong sim, myLong nbSim)
{
    double timeStep = (end_time - start_time) / nb_steps;
    Eigen::VectorXd currentValues = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(spots.data(), spots.size());

    // Generate VDC sequences
    if (vanDerCorputMatrix.rows() != nbSim)
    {
        vanDerCorputMatrix.resize(nbSim, dim);
        for (size_t i = 0; i < nbSim; i++)
        {
            std::vector<double> VDCsequence = gen->GenerateVectorVDC(dim, i);
            for (size_t j = 0; j < dim; j++)
                vanDerCorputMatrix(i, j) = VDCsequence[j];
        }
    }

    // Initialize Trajectories
    for (int i = 0; i < dim; i++)
    {
        trajectory[i] = new SingleTrajectory(start_time, end_time, nb_steps);
        trajectory[i]->AddValue(currentValues[i]);
    }

    for (int i = 0; i < nb_steps; i++)
    {
        // BM Component
        Eigen::VectorXd dW = vanDerCorputMatrix.row(sim) * sqrt(timeStep);

        // Drift Component
        Eigen::VectorXd mapped_drift = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(rates.data(), rates.size());
        Eigen::MatrixXd choleskySquare = choleskyDecomposition * choleskyDecomposition.transpose();
        Eigen::VectorXd milsteinCorrection = 0.5 * choleskySquare.diagonal();

        // Diffusion
        Eigen::VectorXd Z = choleskyDecomposition * dW;
        Eigen::VectorXd X = (mapped_drift - milsteinCorrection) * timeStep + Z + 0.5 * (Z.cwiseProduct(Z));
        currentValues += currentValues.cwiseProduct(X);
        for (int j = 0; j < dim; j++)
        {
            trajectory[j]->AddValue(currentValues[j]);
        }
    }
}
