#pragma once
#include "UniformGenerator.h"

class QuasiGenerator : public UniformGenerator
{
public:
	QuasiGenerator(myLong _current_n = 1);

protected:
	myLong current_n;
};

