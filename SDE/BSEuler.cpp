#include "pch.h"
#include "BSEuler.h"
#include <cmath> 
#include <algorithm>


BSEuler::BSEuler() { }

BSEuler::BSEuler(RandomGenerator* generator, std::vector<double> spots, std::vector<double> rates, Eigen::MatrixXd corrMatrix, int dimension)
    : BlackScholes(generator, spots, rates, corrMatrix, dimension)
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

void BSEuler::Simulate(double start_time, double end_time, size_t nb_steps)
{
    double timeStep = (end_time - start_time) / nb_steps;
    Eigen::VectorXd currentValues = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(spots.data(), spots.size());

    // Initialize  trajectories
    for (int i = 0; i < dim; i++)
    {
        trajectory[i] = new SingleTrajectory(start_time, end_time, nb_steps);
        trajectory[i]->AddValue(currentValues[i]);
    }

    for (int step = 0; step < nb_steps; step++)
    {
        // BM Component
        std::vector<double> dW = gen->GenerateVector(dim) * sqrt(timeStep);
        Eigen::VectorXd mapped_dW = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(dW.data(), dW.size());

        // Drift Component
        std::vector<double> drift = rates * timeStep;
        Eigen::VectorXd mapped_drift = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(drift.data(), drift.size());

        // Generated Trajectories
        Eigen::VectorXd Z = mapped_drift + (choleskyDecomposition * mapped_dW);
        currentValues += currentValues.cwiseProduct(Z);
        for (int i = 0; i < dim; i++)
            trajectory[i]->AddValue(currentValues[i]);
    }
}

void BSEuler::SimulateAntithetic(double start_time, double end_time, size_t nb_steps)
{
    double timeStep = (end_time - start_time) / nb_steps;
    Eigen::VectorXd currentValues = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(spots.data(), spots.size());
    Eigen::VectorXd currentValuesAnti = currentValues;

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
        std::vector<double> dW = gen->GenerateVector(dim) * sqrt(timeStep);
        Eigen::VectorXd mapped_dW = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(dW.data(), dW.size());
        Eigen::VectorXd mapped_dW_anti = -mapped_dW;

        // Drift Component
        std::vector<double> drift = rates * timeStep;
        Eigen::VectorXd mapped_drift = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(drift.data(), drift.size());

        // Generated Trajectoriees
        Eigen::VectorXd Z = mapped_drift + (choleskyDecomposition * mapped_dW);
        Eigen::VectorXd Z_anti = mapped_drift + (choleskyDecomposition * mapped_dW_anti);
        currentValues += currentValues.cwiseProduct(Z);
        currentValuesAnti += currentValuesAnti.cwiseProduct(Z_anti);
        for (int j = 0; j < dim; j++)
        {
            trajectory[j]->AddValue(currentValues[j]);
            trajectoryAntithetic[j]->AddValue(currentValuesAnti[j]);
        }
    }
}

void BSEuler::SimulateQuasiMC(double start_time, double end_time, size_t nb_steps, myLong sim, myLong nbSim)
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
        Eigen::VectorXd dW = vanDerCorputMatrix.row(sim);

        // Drift Component
        std::vector<double> drift = rates * timeStep;
        Eigen::VectorXd mapped_drift = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(drift.data(), drift.size());

        // Generated Trajectoriees
        Eigen::VectorXd Z = mapped_drift + (choleskyDecomposition * dW);
        currentValues += currentValues.cwiseProduct(Z);
        for (int j = 0; j < dim; j++)
        {
            trajectory[j]->AddValue(currentValues[j]);
        }
    }
}

// Helper Function
std::vector<double> operator*(std::vector<double> lhs, double rhs)
{
    size_t dim = lhs.size();
    std::vector<double> result(dim);
    for (size_t i = 0; i < dim; i++)
        result[i] = lhs[i] * rhs;
    return result;
}