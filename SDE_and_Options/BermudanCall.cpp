#include "pch.h"
#include "BermudanCall.h"
#include <eigen-3.4.0/Eigen/Dense>

BermudanCall::BermudanCall()
{
}

BermudanCall::BermudanCall(StochasticProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _exeDates, int _L)
	: Option(_process, _K, _r, _T), exeDates(_exeDates), L(_L)
{
}

double BermudanCall::ComputePrice(int NbSim, bool antithetic)
{

	int nb_exe_dates = exeDates.size();
	Eigen::MatrixXd S_exe_M(NbSim, nb_exe_dates);
	Eigen::MatrixXd S_exe_M_anti(NbSim, nb_exe_dates);
	Eigen::MatrixXd Tau(NbSim, nb_exe_dates);
	Eigen::MatrixXd Tau_anti(NbSim, nb_exe_dates);

	double tk;
	double somme = 0.;
	double somme_anti = 0.;
	double last_value;
	double price;

	v.clear();
	v.resize(NbSim);

	if (antithetic)
	{
		for (int n = 0; n < NbSim; ++n)
		{
			Tau(n, nb_exe_dates - 1) = T;
			Tau_anti(n, nb_exe_dates - 1) = T;
			process->Simulate_Antithetic(0, T, T * 365);

			for (int k = 0; k < nb_exe_dates; k++)
			{
				S_exe_M(n, k) = process->Get_Value(exeDates[k]);
			}

			for (int k = 0; k < nb_exe_dates; k++)
			{
				S_exe_M_anti(n, k) = process->Get_Value_antithetic(exeDates[k]);
			}
		}

		for (int k = nb_exe_dates - 2; k >= 0; k--)
		{
			tk = exeDates[k];
			Eigen::MatrixXd A(NbSim, L);
			Eigen::MatrixXd B(NbSim, 1);
			Eigen::MatrixXd A_anti(NbSim, L);
			Eigen::MatrixXd B_anti(NbSim, 1);

			for (int n = 0; n < NbSim; n++)
			{
				double Tk1 = Tau(n, k + 1);
				double Tk1_anti = Tau_anti(n, k + 1);
				int j = nb_exe_dates - 1;
				while (Tk1 != exeDates[j])
				{
					j--;
				}

				B(n, 0) = std::exp(-r[0] * (Tau(n, k + 1) - tk)) * std::max(S_exe_M(n, j) - K, 0.);

				j = nb_exe_dates - 1;
				while (Tk1_anti != exeDates[j])
				{
					j--;
				}

				B_anti(n, 0) = std::exp(-r[0] * (Tau_anti(n, k + 1) - tk)) * std::max(S_exe_M_anti(n, j) - K, 0.);


				for (int z = 0; z < L; z++)
				{
					A(n, z) = fp(S_exe_M(n, k), z);
				}

				for (int z = 0; z < L; z++)
				{
					A_anti(n, z) = fp(S_exe_M_anti(n, k), z);
				}

			}

			Eigen::VectorXd a(A.colPivHouseholderQr().solve(B));
			Eigen::VectorXd a_anti(A_anti.colPivHouseholderQr().solve(B_anti));

			for (int n = 0; n < NbSim; n++)
			{
				double payoff_actual = std::max(S_exe_M(n, k) - K, 0.);
				double payoff_actual_anti = std::max(S_exe_M_anti(n, k) - K, 0.);
				double expected_payoff = A.row(n) * a;
				double expected_payoff_anti = A_anti.row(n) * a_anti;

				if (payoff_actual >= expected_payoff)
				{
					Tau(n, k) = tk;
				}
				else
				{
					Tau(n, k) = Tau(n, k + 1);
				}

				if (payoff_actual_anti >= expected_payoff_anti)
				{
					Tau_anti(n, k) = tk;
				}
				else
				{
					Tau_anti(n, k) = Tau_anti(n, k + 1);
				}
			}
		}

		for (int n = 0; n < NbSim; n++)
		{
			double T0 = Tau(n, 0);
			int j = nb_exe_dates - 1;
			while (T0 != exeDates[j])
			{
				j--;
			}
			somme += std::exp(-r[0] * T0) * std::max(S_exe_M(n, j) - K, 0.);

			T0 = Tau_anti(n, 0);
			j = nb_exe_dates - 1;
			while (T0 != exeDates[j])
			{
				j--;
			}
			somme_anti += std::exp(-r[0] * T0) * std::max(S_exe_M_anti(n, j) - K, 0.);

			v[n] = (std::exp(-r[0] * T0) * std::max(S_exe_M(n, j) - K, 0.) + std::exp(-r[0] * T0) * std::max(S_exe_M_anti(n, j) - K, 0.)) / 2;
		}

		price = ((somme + somme_anti) / (2 * NbSim));

	}

	else
	{
		for (int n = 0; n < NbSim; ++n)
		{
			Tau(n, nb_exe_dates - 1) = T;
			process->Simulate(0, T, T * 365);


			for (int k = 0; k < nb_exe_dates; k++)
			{
				S_exe_M(n, k) = process->Get_Value(exeDates[k]);
			}

		}

		for (int k = nb_exe_dates - 2; k >= 0; k--)
		{
			tk = exeDates[k];
			Eigen::MatrixXd A(NbSim, L);
			Eigen::MatrixXd B(NbSim, 1);

			for (int n = 0; n < NbSim; n++)
			{
				double Tk1 = Tau(n, k + 1);
				int j = nb_exe_dates - 1;
				while (Tk1 != exeDates[j])
				{
					j--;
				}

				B(n, 0) = std::exp(-r[0] * (Tau(n, k + 1) - tk)) * std::max(S_exe_M(n, j) - K, 0.);


				for (int z = 0; z < L; z++)
				{
					A(n, z) = fp(S_exe_M(n, k), z);
				}

			}

			Eigen::VectorXd a(A.colPivHouseholderQr().solve(B));

			for (int n = 0; n < NbSim; n++)
			{
				double payoff_actual = std::max(S_exe_M(n, k) - K, 0.);
				double expected_payoff = A.row(n) * a;

				if (payoff_actual >= expected_payoff)
				{
					Tau(n, k) = tk;
				}
				else
				{
					Tau(n, k) = Tau(n, k + 1);
				}
			}

		}

		for (int n = 0; n < NbSim; n++)
		{
			double T0 = Tau(n, 0);
			int j = nb_exe_dates - 1;
			while (T0 != exeDates[j])
			{
				j--;
			}
			somme += std::exp(-r[0] * T0) * std::max(S_exe_M(n, j) - K, 0.);
			v[n] = std::exp(-r[0] * T0) * std::max(S_exe_M(n, j) - K, 0.);
		}

		price = (somme / NbSim);
	}

	return price;


}

double BermudanCall::ComputePrice_ControlVariate(int NbSim)
{

	int nb_exe_dates = exeDates.size();
	Eigen::MatrixXd S_exe_M(NbSim, nb_exe_dates);
	Eigen::MatrixXd Tau(NbSim, nb_exe_dates);

	double tk;
	double somme = 0.;
	double last_value;
	double price;

	v.clear();
	v.resize(NbSim);

	for (int n = 0; n < NbSim; ++n)
	{
		Tau(n, nb_exe_dates - 1) = T;
		process->Simulate(0, T, T * 365);

		for (int k = 0; k < nb_exe_dates; k++)
		{
			S_exe_M(n, k) = process->Get_Value(exeDates[k]);
		}

	}
	for (int k = nb_exe_dates - 2; k >= 0; k--)
	{
		tk = exeDates[k];
		Eigen::MatrixXd A(NbSim, L);
		Eigen::MatrixXd B(NbSim, 1);

		for (int n = 0; n < NbSim; n++)
		{
			double Tk1 = Tau(n, k + 1);
			int j = nb_exe_dates - 1;
			while (Tk1 != exeDates[j])
			{
				j--;
			}

			B(n, 0) = std::exp(-r[0] * (Tau(n, k + 1) - tk)) * (std::max(K - S_exe_M(n, j), 0.) + process->Get_Value(0) * std::exp(r[0] * T) - K);


			for (int z = 0; z < L; z++)
			{
				A(n, z) = fp(S_exe_M(n, k), z);
			}

		}

		Eigen::VectorXd a(A.colPivHouseholderQr().solve(B));

		for (int n = 0; n < NbSim; n++)
		{
			double payoff_actual = std::max(K - S_exe_M(n, k), 0.) + process->Get_Value(0) * std::exp(r[0] * T) - K;
			double expected_payoff = A.row(n) * a;

			if (payoff_actual >= expected_payoff)
			{
				Tau(n, k) = tk;
			}
			else
			{
				Tau(n, k) = Tau(n, k + 1);
			}
		}

	}

	for (int n = 0; n < NbSim; n++)
	{
		double T0 = Tau(n, 0);
		int j = nb_exe_dates - 1;
		while (T0 != exeDates[j])
		{
			j--;
		}
		somme += std::exp(-r[0] * T0) * (std::max(K - S_exe_M(n, j), 0.) + process->Get_Value(0) * std::exp(r[0] * T) - K);
		v[n] = std::exp(-r[0] * T0) * (std::max(K - S_exe_M(n, j), 0.) + process->Get_Value(0) * std::exp(r[0] * T) - K);
	}

	price = (somme / NbSim);

	return price;
}

double BermudanCall::ComputePrice_VDC(int NbSim)
{
	int nb_exe_dates = exeDates.size();
	Eigen::MatrixXd S_exe_M(NbSim, nb_exe_dates);
	Eigen::MatrixXd Tau(NbSim, nb_exe_dates);

	double tk;
	double somme = 0.;
	double last_value;
	double price;

	v.clear();
	v.resize(NbSim);

	for (int n = 0; n < NbSim; ++n)
	{
		Tau(n, nb_exe_dates - 1) = T;
		process->Simulate_VDC(0, T, T * 365, n, NbSim);


		for (int k = 0; k < nb_exe_dates; k++)
		{
			S_exe_M(n, k) = process->Get_Value(exeDates[k]);
		}

	}

	for (int k = nb_exe_dates - 2; k >= 0; k--)
	{
		tk = exeDates[k];
		Eigen::MatrixXd A(NbSim, L);
		Eigen::MatrixXd B(NbSim, 1);

		for (int n = 0; n < NbSim; n++)
		{
			double Tk1 = Tau(n, k + 1);
			int j = nb_exe_dates - 1;
			while (Tk1 != exeDates[j])
			{
				j--;
			}

			B(n, 0) = std::exp(-r[0] * (Tau(n, k + 1) - tk)) * std::max(S_exe_M(n, j) - K, 0.);


			for (int z = 0; z < L; z++)
			{
				A(n, z) = fp(S_exe_M(n, k), z);
			}

		}

		Eigen::VectorXd a(A.colPivHouseholderQr().solve(B));

		for (int n = 0; n < NbSim; n++)
		{
			double payoff_actual = std::max(S_exe_M(n, k) - K, 0.);
			double expected_payoff = A.row(n) * a;

			if (payoff_actual >= expected_payoff)
			{
				Tau(n, k) = tk;
			}
			else
			{
				Tau(n, k) = Tau(n, k + 1);
			}
		}

	}

	for (int n = 0; n < NbSim; n++)
	{
		double T0 = Tau(n, 0);
		int j = nb_exe_dates - 1;
		while (T0 != exeDates[j])
		{
			j--;
		}
		somme += std::exp(-r[0] * T0) * std::max(S_exe_M(n, j) - K, 0.);
		v[n] = std::exp(-r[0] * T0) * std::max(S_exe_M(n, j) - K, 0.);
	}

	price = (somme / NbSim);
	return price;
}