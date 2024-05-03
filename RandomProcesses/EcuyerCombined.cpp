#include "pch.h"
#include "EcuyerCombined.h"
#include "iostream"

EcuyerCombined::EcuyerCombined()
	: Generator1(10, 40014, 0, 2147483563), Generator2(15, 40692, 0, 2147483399)
{

}

double EcuyerCombined::Generate()
{
	myLong current_number1;
	myLong current_number2;

	current_number1 = Generator1.Generate() * Generator1.get_Modulus();
	current_number2 = Generator2.Generate() * Generator2.get_Modulus();
	currentNumber = (current_number1 - current_number2) & (Generator1.get_Modulus() - 1);

	if (currentNumber > 0)
	{
		return (double)currentNumber / Generator1.get_Modulus();
	}
	else if (currentNumber = 0)
	{
		return (double)(Generator1.get_Modulus() - 1) / Generator1.get_Modulus();
	}
	else
	{
		return (double)(-1) * currentNumber / Generator1.get_Modulus();
	}
}