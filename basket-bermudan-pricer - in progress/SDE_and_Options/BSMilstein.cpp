#include "pch.h"
#include "BSMilstein.h"
#include <math.h>
#include <algorithm>
#include <random>

using namespace Eigen;


BSMilstein::BSMilstein()
{
}

BSMilstein::BSMilstein(RandomGenerator* _gen, std::vector<double> _s, std::vector<double> _r, Eigen::MatrixXd _VarCovar, int dim)
    : BlackScholes(_gen, _s, _r, _VarCovar, dim)
{
    EigenSolver<MatrixXd> es(VarCovar);
    if (VarCovar.determinant() == 0)
    {
        MatrixXd D = es.pseudoEigenvalueMatrix();
        MatrixXd P = es.pseudoEigenvectors();
        B = P * D;
    }
    else
    {
        B = VarCovar.llt().matrixL();
    }

}

void BSMilstein::Simulate(double start_time, double end_time, size_t nb_steps)
{

    double dt = (end_time - start_time) / nb_steps;
    Eigen::VectorXd last = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(s.data(), s.size());


    for (int i = 0; i < dim; i++)
    {
        trajectory[i] = new SingleTrajectory(start_time, end_time, nb_steps);
        trajectory[i]->AddValue(last[i]);
    }


    for (int i = 0; i < nb_steps; i++)
    {  
        std::vector<double> dW = gen->GenerateVector(dim);
        Eigen::VectorXd dW_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(dW.data(), dW.size())* pow(dt, 0.5);

        Eigen::VectorXd R_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(r.data(), r.size());

        Eigen::MatrixXd B_square = B * B;

        Eigen::VectorXd T = 0.5 * VarCovar.diagonal();

        Eigen::VectorXd Z = B * dW_M;

        Eigen::VectorXd X = (R_M - T) * dt + Z + 0.5 * (Z.cwiseProduct(Z));

        Eigen::VectorXd next = last + last.cwiseProduct(X);
        last = next;

        for (int j = 0; j < dim; j++)
        {
            trajectory[j]->AddValue(last[j]);
        }

    }
}

void BSMilstein::SimulateAntithetic(double start_time, double end_time, size_t nb_steps)
{
    double dt = (end_time - start_time) / nb_steps;
    Eigen::VectorXd last = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(s.data(), s.size());
    Eigen::VectorXd last_anti = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(s.data(), s.size());

    for (int i = 0; i < dim; i++)
    {
        trajectory[i] = new SingleTrajectory(start_time, end_time, nb_steps);
        trajectory[i]->AddValue(last[i]);

        trajectoryAntithetic[i] = new SingleTrajectory(start_time, end_time, nb_steps);
        trajectoryAntithetic[i]->AddValue(last_anti[i]);
    }

    for (int i = 0; i < nb_steps; i++)
    {

        std::vector<double> dW = gen->GenerateVector(dim);
        Eigen::VectorXd dW_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(dW.data(), dW.size())* pow(dt, 0.5);
        Eigen::VectorXd R_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(r.data(), r.size());
        Eigen::VectorXd dW_M_anti = -dW_M;

        Eigen::MatrixXd B_square = B * B;

        Eigen::VectorXd T = 0.5 * VarCovar.diagonal();

        Eigen::VectorXd Z = B * dW_M;
        Eigen::VectorXd Z_anti = B * dW_M_anti;

        Eigen::VectorXd X = (R_M - T) * dt + Z + 0.5 * (Z.cwiseProduct(Z));
        Eigen::VectorXd X_anti = (R_M - T) * dt + Z_anti + 0.5 * (Z_anti.cwiseProduct(Z_anti));

        Eigen::VectorXd next = last + last.cwiseProduct(X);
        Eigen::VectorXd next_anti = last_anti + last_anti.cwiseProduct(X_anti);

        last = next;
        last_anti = next_anti;

        for (int j = 0; j < dim; j++)
        {
            trajectory[j]->AddValue(last[j]);
            trajectoryAntithetic[j]->AddValue(last_anti[j]);
        }

    }
}

void BSMilstein::SimulateQuasiMC(double start_time, double end_time, size_t nb_steps, myLong sim, myLong nbSim)
{
    double dt = (end_time - start_time) / nb_steps;
    Eigen::VectorXd last = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(s.data(), s.size());
    Eigen::MatrixXd x;
    std::vector<double> y;

    if (M_VDC.rows() != nbSim)
    {
        M_VDC.resize(nbSim, dim);
        for (size_t i = 0; i < nbSim; i++)
        {
            y = gen->GenerateVectorVDC(dim, i);
            for (size_t j = 0; j < dim; j++)
            {
                M_VDC(i, j) = y[j];
            }
        }

    }

    for (int i = 0; i < dim; i++)
    {
        trajectory[i] = new SingleTrajectory(start_time, end_time, nb_steps);
        trajectory[i]->AddValue(last[i]);
    }

    Eigen::VectorXd dW_M = M_VDC.row(sim) * pow(dt, 0.5);

    Eigen::VectorXd R_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(r.data(), r.size());

    Eigen::MatrixXd B_square = B * B;

    Eigen::VectorXd T = 0.5 * VarCovar.diagonal();

    Eigen::VectorXd Z = B * dW_M;

    Eigen::VectorXd X = (R_M - T) * dt + Z + 0.5 * (Z.cwiseProduct(Z));

    Eigen::VectorXd next = last + last.cwiseProduct(X);
    last = next;

    for (int j = 0; j < dim; j++)
    {
        trajectory[j]->AddValue(last[j]);
    }

    for (int i = 0; i < nb_steps; i++)
    {
        x = M_VDC;

        for (unsigned int j = 0; j < dim; j++)
        {
            auto rng = std::default_random_engine{ i * (j + 1) };
            std::shuffle(x.col(j).begin(), x.col(j).end(), rng); //Needed permutation for quasi-monte carlo
        }
        Eigen::VectorXd dW_M = x.row(sim) * pow(dt, 0.5);

        Eigen::VectorXd R_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(r.data(), r.size());

        Eigen::MatrixXd B_square = B * B;

        Eigen::VectorXd T = 0.5 * VarCovar.diagonal();

        Eigen::VectorXd Z = B * dW_M;

        Eigen::VectorXd X = (R_M - T) * dt + Z + 0.5 * (Z.cwiseProduct(Z));

        Eigen::VectorXd next = last + last.cwiseProduct(X);
        last = next;

        for (int j = 0; j < dim; j++)
        {
            trajectory[j]->AddValue(last[j]);
        }

    }
}
