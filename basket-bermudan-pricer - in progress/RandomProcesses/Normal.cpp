#include "pch.h"
#include "Normal.h"
#include <stdexcept>
#include <cmath>

const long double PI = 3.141592653589793238L;

/*-- Base Class --*/
Normal::Normal() { }

Normal::Normal(double _mean, double _var, UniformGenerator* _gen) : mean(_mean), var(_var), ContinuousGenerator(_gen) { }

/*-- Normal Box Muller--*/
NormalBoxMuller::NormalBoxMuller() { }

NormalBoxMuller::NormalBoxMuller(double _mean, double _var, UniformGenerator* _gen) : Normal(_mean, _var, _gen), SecondNormal(0) { }

double NormalBoxMuller::Generate()
{
	if (NewSimulation)
	{
		double firstUniform = generator->Generate();
		double secondUniform = generator->Generate();

		double R = pow(-2 * log(firstUniform), 0.5);
		double theta = 2 * PI * secondUniform;
		double firstNormal = R * cos(theta);

		SecondNormal = R * sin(theta);
		NewSimulation = false;

		return firstNormal * pow(var, 0.5) + mean;
	}
	else
	{
		NewSimulation = true;
		return SecondNormal * pow(var, 0.5) + mean;
	}
}

/*-- Normal CLT--*/
NormalCLT::NormalCLT() { }

NormalCLT::NormalCLT(double _mean, double _var, UniformGenerator* _gen) : Normal(_mean, _var, _gen) { }

double NormalCLT::Generate()
{
	double sumUniform = 0;
	for (int i = 0; i < 12; i++)
		sumUniform += generator->Generate();

	return (sumUniform - 6) * pow(var, 0.5) + mean;
}

/*-- Normal InvCDF--*/
NormalInvCdf::NormalInvCdf() { }

NormalInvCdf::NormalInvCdf(double _mean, double _var, UniformGenerator* _gen) : Normal(_mean, _var, _gen) { }

double RationalApproximation(double t)
{
    // Abramowitz and Stegun formula
    double c[] = { 2.515517, 0.802853, 0.010328 };
    double d[] = { 1.432788, 0.189269, 0.001308 };
    return t - ((c[2] * t + c[1]) * t + c[0]) /
        (((d[2] * t + d[1]) * t + d[0]) * t + 1.0);
}

double NormalCDFInverse(double p)
{
    if (p <= 0.0 || p >= 1.0)
        throw std::invalid_argument("p must be between 0 and 1.");
    if (p < 0.5)
        return -RationalApproximation(sqrt(-2.0 * log(p)));
    else
        return RationalApproximation(sqrt(-2.0 * log(1 - p)));
}

double NormalInvCdf::Generate()
{
    double Uniform = generator->Generate();
    return NormalCDFInverse(Uniform);
}