#include "pch.h"
#include "BermudanBasketCall.h"
#include "ClosedFormulas.h"
#include <numeric>
#include <iterator>

Bermudan::Bermudan()
{
}

Bermudan::Bermudan(StochasticProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights, std::vector<double> _exeDates, int _L)
	: Option(_process, _K, _r, _T), exeDates(_exeDates), L(_L), weights(_weights)
{
}

Bermudan::Bermudan(StochasticProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights, std::vector<double> _exeDates, std::vector<double> _S, Eigen::MatrixXd _VarCovar, int _L)
	: Option(_process, _K, _r, _T), exeDates(_exeDates), L(_L), weights(_weights), S(_S), VarCovar(_VarCovar)
{
}

double Bermudan::ComputePrice(int NbSim, bool antithetic)
{

	int nb_exe_dates = exeDates.size();
	Eigen::MatrixXd S_exe_M(NbSim, nb_exe_dates);
	Eigen::MatrixXd S_exe_M_anti(NbSim, nb_exe_dates);
	Eigen::MatrixXd Tau(NbSim, nb_exe_dates);
	Eigen::MatrixXd Tau_anti(NbSim, nb_exe_dates);
	Eigen::VectorXd weights_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

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
				std::vector<double> Snk = process->Get_ValueND(exeDates[k]);
				Eigen::VectorXd Snk_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Snk.data(), Snk.size());
				S_exe_M(n, k) = weights_V.transpose() * Snk_V;
			}

			for (int k = 0; k < nb_exe_dates; k++)
			{
				std::vector<double> Snk_anti = process->Get_ValueND_antithetic(exeDates[k]);
				Eigen::VectorXd Snk_V_anti = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Snk_anti.data(), Snk_anti.size());
				S_exe_M_anti(n, k) = weights_V.transpose() * Snk_V_anti;
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
				double expected_payoff = A.row(n) * a;

				double payoff_actual_anti = std::max(S_exe_M_anti(n, k) - K, 0.);
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
				std::vector<double> Snk = process->Get_ValueND(exeDates[k]);
				Eigen::VectorXd Snk_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Snk.data(), Snk.size());
				S_exe_M(n, k) = weights_V.transpose() * Snk_V;
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

double Bermudan::ComputePrice_ControlVariate(int NbSim)
{

	int nb_exe_dates = exeDates.size();
	Eigen::MatrixXd S_exe_M(NbSim, nb_exe_dates);
	Eigen::MatrixXd S_exe_M_L(NbSim, nb_exe_dates);
	Eigen::MatrixXd Tau(NbSim, nb_exe_dates);
	Eigen::VectorXd weights_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

	double tk;
	double somme = 0.;
	double last_value;
	double price;

	v.clear();
	v.resize(NbSim);
	double E_Y = Compute_E_Ybis(S, weights, K, r[0], VarCovar, T);


	for (int n = 0; n < NbSim; ++n)
	{
		Tau(n, nb_exe_dates - 1) = T;
		process->Simulate(0, T, T * 365);

		for (int k = 0; k < nb_exe_dates; k++)
		{
			std::vector<double> Snk = process->Get_ValueND(exeDates[k]);
			std::vector<double> Snk_L(Snk.size());
			for (int j = 0; j < S.size(); j++)
			{
				Snk_L[j] = log(Snk[j]);
			}
			Eigen::VectorXd Snk_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Snk.data(), Snk.size());
			Eigen::VectorXd Snk_V_L = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Snk_L.data(), Snk_L.size());
			S_exe_M(n, k) = weights_V.transpose() * Snk_V;
			S_exe_M_L(n, k) = weights_V.transpose() * Snk_V_L;
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

			B(n, 0) = std::exp(-r[0] * (Tau(n, k + 1) - tk)) * (std::max(S_exe_M(n, j) - K, 0.) - std::max(exp(S_exe_M_L(n, j)) - K, 0.) + E_Y);

			for (int z = 0; z < L; z++)
			{
				A(n, z) = fp(S_exe_M(n, k), z);
			}

		}

		Eigen::VectorXd a(A.colPivHouseholderQr().solve(B));

		for (int n = 0; n < NbSim; n++)
		{
			double payoff_actual = std::max(S_exe_M(n, k) - K, 0.) - std::max(exp(S_exe_M_L(n, k)) - K, 0.) + E_Y;
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
		somme += std::exp(-r[0] * T0) * (std::max(S_exe_M(n, j) - K, 0.) - std::max(exp(S_exe_M_L(n, j)) - K, 0.) + E_Y);
		v[n] = std::exp(-r[0] * T0) * (std::max(S_exe_M(n, j) - K, 0.) - std::max(exp(S_exe_M_L(n, j)) - K, 0.) + E_Y);
	}

	price = (somme / NbSim);
	return price;
}


double Compute_E_Ybis(std::vector<double> S, std::vector<double> weights, double K, double r, Eigen::MatrixXd VarCovar, double T)
{

	double S_Y = 1;
	double R_Y;
	double vol_Y;

	Eigen::VectorXd weights_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());
	Eigen::MatrixXd B;

	Eigen::EigenSolver<Eigen::MatrixXd> es(VarCovar);

	if (VarCovar.determinant() == 0)
	{
		Eigen::MatrixXd D = es.pseudoEigenvalueMatrix();
		Eigen::MatrixXd P = es.pseudoEigenvectors();
		B = P * D;
	}
	else
	{
		B = VarCovar.llt().matrixL();
	}

	for (int i = 0; i < S.size(); i++)
	{
		S_Y *= pow(S[i], weights[i]);
	}


	Eigen::MatrixXd sigma2 = B * B;
	Eigen::VectorXd sigmai = sigma2.colwise().sum();
	double R1 = 0.5 * weights_M.transpose() * sigmai;
	double R2 = 0.5 * weights_M.transpose() * B * B.transpose() * weights_M;
	R_Y = r - R1 + R2;


	vol_Y = pow(weights_M.transpose() * B * B.transpose() * weights_M, 0.5);

	double price = bs_price_call(S_Y, K, vol_Y, T, R_Y);

	return price;

}

double Bermudan::ComputePrice_VDC(int NbSim)
{
	int nb_exe_dates = exeDates.size();
	Eigen::MatrixXd S_exe_M(NbSim, nb_exe_dates);
	Eigen::MatrixXd Tau(NbSim, nb_exe_dates);
	Eigen::VectorXd weights_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

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
			std::vector<double> Snk = process->Get_ValueND(exeDates[k]);
			Eigen::VectorXd Snk_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Snk.data(), Snk.size());
			S_exe_M(n, k) = weights_V.transpose() * Snk_V;
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