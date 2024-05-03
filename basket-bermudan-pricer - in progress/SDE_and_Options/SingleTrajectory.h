#pragma once
#include <vector>

class SingleTrajectory
{
public:
	SingleTrajectory();
	SingleTrajectory(double _start, double _end, size_t _steps);
	void AddValue(double val);
	const double GetValue(double time);

protected:
	double step;
	double start;
	double end;
	size_t steps;

	std::vector<double> Values;
	std::vector<double> Times;
};