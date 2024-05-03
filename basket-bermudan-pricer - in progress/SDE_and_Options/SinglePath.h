#pragma once
#include <vector>
#include <iostream>

//Store random process path 

class SinglePath
{
public:
	SinglePath();
	SinglePath(double _start, double _end, size_t _nbSteps);
	void AddValue(double val);
	const double GetValue(double time);
	void print_vector();

protected:
	double timeStep;
	std::vector<double> Values;
	std::vector<double> Times;
	double start;
	double end;
	size_t nbSteps;
};

