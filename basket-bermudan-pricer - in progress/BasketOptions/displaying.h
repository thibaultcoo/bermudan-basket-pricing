#pragma once

#include <iostream>
#include "pricing.h"

class Display 
{
public:
    Display(Priceable* pric);
    void display(Priceable* pric, bool quasiMC);
};