#pragma once
#include "UniformGenerator.h"

class QuasiGenerator : public UniformGenerator
{
	protected:
		myLong current_n;

	public:
		QuasiGenerator(myLong _current_n = 1);
};