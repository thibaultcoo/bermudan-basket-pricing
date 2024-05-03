#include "pch.h"
#include "LinearCongruential.h"
#include "iostream"

LinearCongruential::LinearCongruential()
	: PseudoGenerator()
{
}


LinearCongruential::LinearCongruential(myLong _seed, myLong _multiplier, myLong _increment, myLong _modulus)
	: PseudoGenerator(_seed), Multiplier(_multiplier), Increment(_increment), Modulus(_modulus)
{

}

double LinearCongruential::Generate()
{
	currentNumber = (Multiplier * currentNumber + Increment) % Modulus;
	return (double)currentNumber / Modulus;
}

myLong LinearCongruential::get_Modulus()
{
	return Modulus;
}