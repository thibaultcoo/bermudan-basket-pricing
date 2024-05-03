#include "pch.h"
#include "BSEuler.h"
#include <math.h> 
#include <algorithm>
#include <random>

using namespace Eigen;


BSEuler::BSEuler()
{
}

BSEuler::BSEuler(RandomGenerator* _gen, std::vector<double> _s, std::vector<double> _r, Eigen::MatrixXd _VarCovar, int dim)
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

void BSEuler::Simulate(double start_time, double end_time, size_t nb_steps)
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

        std::vector<double> dW = gen->GenerateVector(dim) * pow(dt, 0.5);
        Eigen::VectorXd dW_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(dW.data(), dW.size());

        std::vector<double> M = r * dt;

        Eigen::VectorXd M_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(M.data(), M.size());
        Eigen::VectorXd Z = M_M + (B * dW_M);
        Eigen::VectorXd next = last + last.cwiseProduct(Z);
        last = next;
        for (int j = 0; j < dim; j++)
        {
            trajectory[j]->AddValue(last[j]);
        }

    }
}

void BSEuler::SimulateAntithetic(double start_time, double end_time, size_t nb_steps)
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

        std::vector<double> dW = gen->GenerateVector(dim) * pow(dt, 0.5);
        Eigen::VectorXd dW_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(dW.data(), dW.size());
        Eigen::VectorXd dW_M_anti = -dW_M;
        std::vector<double> M = r * dt;

        Eigen::VectorXd M_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(M.data(), M.size());
        Eigen::VectorXd Z = M_M + (B * dW_M);
        Eigen::VectorXd Z_anti = M_M + (B * dW_M_anti);

        Eigen::VectorXd next = last + last.cwiseProduct(Z);
        last = next;

        Eigen::VectorXd next_anti = last_anti + last_anti.cwiseProduct(Z_anti);
        last_anti = next_anti;

        for (int j = 0; j < dim; j++)
        {
            trajectory[j]->AddValue(last[j]);
            trajectoryAntithetic[j]->AddValue(last_anti[j]);
        }

    }
}

void BSEuler::SimulateQuasiMC(double start_time, double end_time, size_t nb_steps, myLong sim, myLong nbSim)
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
    std::vector<double> M = r * dt;

    Eigen::VectorXd M_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(M.data(), M.size());
    Eigen::VectorXd Z = M_M + (B * dW_M);
    Eigen::VectorXd next = last + last.cwiseProduct(Z);
    last = next;

    for (int j = 0; j < dim; j++)
    {
        trajectory[j]->AddValue(last[j]);
    }


    for (int i = 1; i < nb_steps; i++)
    {
        x = M_VDC;
        for (unsigned int j = 0; j < dim; j++)
        {
            auto rng = std::default_random_engine{ i * (j + 1) };
            std::shuffle(x.col(j).begin(), x.col(j).end(), rng); //Needed permutation for quasi-monte carlo
        }

        Eigen::VectorXd dW = x.row(sim) * pow(dt, 0.5);
        std::vector<double> M = r * dt;
        Eigen::VectorXd M_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(M.data(), M.size());
        Eigen::VectorXd Z = M_M + (B * dW);
        Eigen::VectorXd next = last + last.cwiseProduct(Z);
        last = next;
        for (int j = 0; j < dim; j++)
        {
            trajectory[j]->AddValue(last[j]);
        }

    }
}

std::vector<double> operator*(std::vector<double> lhs, double rhs)
{
    size_t dim = lhs.size();
    std::vector<double> result(dim);
    for (size_t i = 0; i < dim; i++)
    {
        result[i] = lhs[i] * rhs;
    }
    return result;
}