#include "pch.h"
#include "VanDerCorput.h"
#include <vector>

VanDerCorput::VanDerCorput(int _base, myLong _currentNumber) : QuasiGenerator(_currentNumber), base(_base) { }

double VanDerCorput::Generate()
{
    std::vector<int> b = InversePAdicExpansion(current_n, base);
    double phi = 0;
    for (int k = 0; k < b.size(); k++)
    {
        phi += b[k] / pow(base, (k + 1));
    }
    current_n = current_n + 1;
    return phi;
}

std::vector<int> InversePAdicExpansion(int n, int base)
{
    std::vector<int> b;

    while (n > 0)
    {
        int d = n % base;
        b.push_back(d);
        n = n / base;
    }

    return b;
}
