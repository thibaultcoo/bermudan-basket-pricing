#include "BermudanBasket.h"

BermudanBasket::BermudanBasket() { }

BermudanBasket::BermudanBasket(StochasticProcess* process, double strike, std::vector<double> rates, double maturity, std::vector<double> weights, std::vector<double> exerciseDates, int L)
	: Options(process, strike, rates, maturity), exerciseDates(exerciseDates), L(L), weights(weights) { }

BermudanBasket::BermudanBasket(StochasticProcess* process, double strike, std::vector<double> rates, double maturity, std::vector<double> weights, std::vector<double> exerciseDates, std::vector<double> spots, Eigen::MatrixXd corrMatrix, int L)
	: Options(process, strike, rates, maturity), exerciseDates(exerciseDates), L(L), weights(weights), spots(spots), corrMatrix(corrMatrix) { }

double BermudanBasket::price(int nbSim)
{
	int nb_exerciseDates = exerciseDates.size();
	Eigen::MatrixXd S_exe_M(nbSim, nb_exerciseDates);
	Eigen::MatrixXd Tau(nbSim, nb_exerciseDates);
	Eigen::VectorXd weights_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

	double t_k;
	double somme = 0.;
	double price;

	paths_prices.clear();
	paths_prices.resize(nbSim);

	for (int n = 0; n < nbSim; ++n)
	{
		Tau(n, nb_exerciseDates - 1) = maturity;
		process->Simulate(0, maturity, maturity * 365);

		for (int k = 0; k < nb_exerciseDates; k++)
		{
			std::vector<double> Snk = process->Get_ValueND(exerciseDates[k]);
			Eigen::VectorXd Snk_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Snk.data(), Snk.size());
			S_exe_M(n, k) = weights_V.transpose() * Snk_V;
		}
	}

	for (int k = nb_exerciseDates - 2; k >= 0; k--)
	{
		t_k = exerciseDates[k];
		Eigen::MatrixXd A(nbSim, L);
		Eigen::MatrixXd B(nbSim, 1);

		for (int n = 0; n < nbSim; n++)
		{
			double Tk1 = Tau(n, k + 1);
			int j = nb_exerciseDates - 1;
			while (Tk1 != exerciseDates[j])
				j--;

			B(n, 0) = std::exp(-rates[0] * (Tau(n, k + 1) - t_k)) * std::max(S_exe_M(n, j) - strike, 0.);


			for (int z = 0; z < L; z++)
				A(n, z) = fp(S_exe_M(n, k), z);
		}

		Eigen::VectorXd a(A.colPivHouseholderQr().solve(B));

		for (int n = 0; n < nbSim; n++)
		{
			double payoff_actual = std::max(S_exe_M(n, k) - strike, 0.);
			double expected_payoff = A.row(n) * a;

			if (payoff_actual >= expected_payoff)
				Tau(n, k) = t_k;
			else
				Tau(n, k) = Tau(n, k + 1);
		}

	}

	for (int n = 0; n < nbSim; n++)
	{
		double T0 = Tau(n, 0);
		int j = nb_exerciseDates - 1;
		while (T0 != exerciseDates[j])
			j--;

		somme += std::exp(-rates[0] * T0) * std::max(S_exe_M(n, j) - strike, 0.);
		paths_prices[n] = std::exp(-rates[0] * T0) * std::max(S_exe_M(n, j) - strike, 0.);
	}

	price = (somme / nbSim);

	return price;
}

double BermudanBasket::priceAntithetic(int nbSim)
{
	int nb_exerciseDates = exerciseDates.size();
	Eigen::MatrixXd S_exe_M(nbSim, nb_exerciseDates);
	Eigen::MatrixXd S_exe_M_anti(nbSim, nb_exerciseDates);
	Eigen::MatrixXd Tau(nbSim, nb_exerciseDates);
	Eigen::MatrixXd Tau_anti(nbSim, nb_exerciseDates);
	Eigen::VectorXd weights_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

	double t_k;
	double somme = 0.;
	double somme_anti = 0.;
	double last_value;

	paths_prices.clear();
	paths_prices.resize(nbSim);

	for (int n = 0; n < nbSim; ++n)
	{
		Tau(n, nb_exerciseDates - 1) = maturity;
		Tau_anti(n, nb_exerciseDates - 1) = maturity;
		process->Simulate_Antithetic(0, maturity, maturity * 365);

		for (int k = 0; k < nb_exerciseDates; k++)
		{
			std::vector<double> Snk = process->Get_ValueND(exerciseDates[k]);
			Eigen::VectorXd Snk_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Snk.data(), Snk.size());
			S_exe_M(n, k) = weights_V.transpose() * Snk_V;
		}

		for (int k = 0; k < nb_exerciseDates; k++)
		{
			std::vector<double> Snk_anti = process->Get_ValueND_antithetic(exerciseDates[k]);
			Eigen::VectorXd Snk_V_anti = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Snk_anti.data(), Snk_anti.size());
			S_exe_M_anti(n, k) = weights_V.transpose() * Snk_V_anti;
		}

	}

	for (int k = nb_exerciseDates - 2; k >= 0; k--)
	{
		t_k = exerciseDates[k];
		Eigen::MatrixXd A(nbSim, L);
		Eigen::MatrixXd B(nbSim, 1);
		Eigen::MatrixXd A_anti(nbSim, L);
		Eigen::MatrixXd B_anti(nbSim, 1);

		for (int n = 0; n < nbSim; n++)
		{
			double Tk1 = Tau(n, k + 1);
			double Tk1_anti = Tau_anti(n, k + 1);
			int j = nb_exerciseDates - 1;
			while (Tk1 != exerciseDates[j])
				j--;

			B(n, 0) = std::exp(-rates[0] * (Tau(n, k + 1) - t_k)) * std::max(S_exe_M(n, j) - strike, 0.);

			j = nb_exerciseDates - 1;
			while (Tk1_anti != exerciseDates[j])
				j--;

			B_anti(n, 0) = std::exp(-rates[0] * (Tau_anti(n, k + 1) - t_k)) * std::max(S_exe_M_anti(n, j) - strike, 0.);

			for (int z = 0; z < L; z++)
				A(n, z) = fp(S_exe_M(n, k), z);

			for (int z = 0; z < L; z++)
				A_anti(n, z) = fp(S_exe_M_anti(n, k), z);

		}

		Eigen::VectorXd a(A.colPivHouseholderQr().solve(B));
		Eigen::VectorXd a_anti(A_anti.colPivHouseholderQr().solve(B_anti));

		for (int n = 0; n < nbSim; n++)
		{
			double payoff_actual = std::max(S_exe_M(n, k) - strike, 0.);
			double expected_payoff = A.row(n) * a;

			double payoff_actual_anti = std::max(S_exe_M_anti(n, k) - strike, 0.);
			double expected_payoff_anti = A_anti.row(n) * a_anti;

			if (payoff_actual >= expected_payoff)
				Tau(n, k) = t_k;
			else
				Tau(n, k) = Tau(n, k + 1);

			if (payoff_actual_anti >= expected_payoff_anti)
				Tau_anti(n, k) = t_k;
			else
				Tau_anti(n, k) = Tau_anti(n, k + 1);
		}
	}

	for (int n = 0; n < nbSim; n++)
	{
		double T0 = Tau(n, 0);
		int j = nb_exerciseDates - 1;
		while (T0 != exerciseDates[j])
			j--;
		somme += std::exp(-rates[0] * T0) * std::max(S_exe_M(n, j) - strike, 0.);

		T0 = Tau_anti(n, 0);
		j = nb_exerciseDates - 1;
		while (T0 != exerciseDates[j])
			j--;

		somme_anti += std::exp(-rates[0] * T0) * std::max(S_exe_M_anti(n, j) - strike, 0.);
		paths_prices[n] = (std::exp(-rates[0] * T0) * std::max(S_exe_M(n, j) - strike, 0.) + std::exp(-rates[0] * T0) * std::max(S_exe_M_anti(n, j) - strike, 0.)) / 2;
	}

	double price = ((somme + somme_anti) / (2 * nbSim));

	return price;
}

double BermudanBasket::priceControlVariate(int nbSim)
{

	int nb_exerciseDates = exerciseDates.size();
	Eigen::MatrixXd S_exe_M(nbSim, nb_exerciseDates);
	Eigen::MatrixXd S_exe_M_L(nbSim, nb_exerciseDates);
	Eigen::MatrixXd Tau(nbSim, nb_exerciseDates);
	Eigen::VectorXd weights_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

	double tk;
	double somme = 0.;
	double last_value;

	paths_prices.clear();
	paths_prices.resize(nbSim);
	double E_Y = Compute_E_Ybis(spots, weights, strike, rates[0], corrMatrix, maturity);


	for (int n = 0; n < nbSim; ++n)
	{
		Tau(n, nb_exerciseDates - 1) = maturity;
		process->Simulate(0, maturity, maturity * 365);

		for (int k = 0; k < nb_exerciseDates; k++)
		{
			std::vector<double> Snk = process->Get_ValueND(exerciseDates[k]);
			std::vector<double> Snk_L(Snk.size());
			for (int j = 0; j < spots.size(); j++)
				Snk_L[j] = log(Snk[j]);

			Eigen::VectorXd Snk_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Snk.data(), Snk.size());
			Eigen::VectorXd Snk_V_L = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Snk_L.data(), Snk_L.size());
			S_exe_M(n, k) = weights_V.transpose() * Snk_V;
			S_exe_M_L(n, k) = weights_V.transpose() * Snk_V_L;
		}
	}


	for (int k = nb_exerciseDates - 2; k >= 0; k--)
	{
		tk = exerciseDates[k];
		Eigen::MatrixXd A(nbSim, L);
		Eigen::MatrixXd B(nbSim, 1);

		for (int n = 0; n < nbSim; n++)
		{
			double Tk1 = Tau(n, k + 1);
			int j = nb_exerciseDates - 1;
			while (Tk1 != exerciseDates[j])
				j--;

			B(n, 0) = std::exp(-rates[0] * (Tau(n, k + 1) - tk)) * (std::max(S_exe_M(n, j) - strike, 0.) - std::max(exp(S_exe_M_L(n, j)) - strike, 0.) + E_Y);

			for (int z = 0; z < L; z++)
				A(n, z) = fp(S_exe_M(n, k), z);
		}

		Eigen::VectorXd a(A.colPivHouseholderQr().solve(B));

		for (int n = 0; n < nbSim; n++)
		{
			double payoff_actual = std::max(S_exe_M(n, k) - strike, 0.) - std::max(exp(S_exe_M_L(n, k)) - strike, 0.) + E_Y;
			double expected_payoff = A.row(n) * a;

			if (payoff_actual >= expected_payoff)
				Tau(n, k) = tk;
			else
				Tau(n, k) = Tau(n, k + 1);
		}
	}

	for (int n = 0; n < nbSim; n++)
	{
		double T0 = Tau(n, 0);
		int j = nb_exerciseDates - 1;
		while (T0 != exerciseDates[j])
			j--;

		somme += std::exp(-rates[0] * T0) * (std::max(S_exe_M(n, j) - strike, 0.) - std::max(exp(S_exe_M_L(n, j)) - strike, 0.) + E_Y);
		paths_prices[n] = std::exp(-rates[0] * T0) * (std::max(S_exe_M(n, j) - strike, 0.) - std::max(exp(S_exe_M_L(n, j)) - strike, 0.) + E_Y);
	}

	double price = (somme / nbSim);
	
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

	double price = BS_Call(S_Y, K, vol_Y, T, R_Y);

	return price;
}

double BermudanBasket::priceVDC(int nbSim)
{
	int nb_exerciseDates = exerciseDates.size();
	Eigen::MatrixXd S_exe_M(nbSim, nb_exerciseDates);
	Eigen::MatrixXd Tau(nbSim, nb_exerciseDates);
	Eigen::VectorXd weights_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

	double tk;
	double somme = 0.;
	double last_value;

	paths_prices.clear();
	paths_prices.resize(nbSim);

	for (int n = 0; n < nbSim; ++n)
	{
		Tau(n, nb_exerciseDates - 1) = maturity;
		process->Simulate_VDC(0, maturity, maturity * 365, n, nbSim);

		for (int k = 0; k < nb_exerciseDates; k++)
		{
			std::vector<double> Snk = process->Get_ValueND(exerciseDates[k]);
			Eigen::VectorXd Snk_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Snk.data(), Snk.size());
			S_exe_M(n, k) = weights_V.transpose() * Snk_V;
		}
	}

	for (int k = nb_exerciseDates - 2; k >= 0; k--)
	{
		tk = exerciseDates[k];
		Eigen::MatrixXd A(nbSim, L);
		Eigen::MatrixXd B(nbSim, 1);

		for (int n = 0; n < nbSim; n++)
		{
			double Tk1 = Tau(n, k + 1);
			int j = nb_exerciseDates - 1;
			while (Tk1 != exerciseDates[j])
				j--;

			B(n, 0) = std::exp(-rates[0] * (Tau(n, k + 1) - tk)) * std::max(S_exe_M(n, j) - strike, 0.);

			for (int z = 0; z < L; z++)
				A(n, z) = fp(S_exe_M(n, k), z);
		}

		Eigen::VectorXd a(A.colPivHouseholderQr().solve(B));

		for (int n = 0; n < nbSim; n++)
		{
			double payoff_actual = std::max(S_exe_M(n, k) - strike, 0.);
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
	for (int n = 0; n < nbSim; n++)
	{
		double T0 = Tau(n, 0);
		int j = nb_exerciseDates - 1;
		while (T0 != exerciseDates[j])
			j--;

		somme += std::exp(-rates[0] * T0) * std::max(S_exe_M(n, j) - strike, 0.);
		paths_prices[n] = std::exp(-rates[0] * T0) * std::max(S_exe_M(n, j) - strike, 0.);
	}

	double price = (somme / nbSim);
	
	return price;
}