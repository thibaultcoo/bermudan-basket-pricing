#pragma once

#include <iostream>
#include "pricing.h"

class Display
{
public:
    // Constructor to initialize with a Priceable object
    explicit Display(Priceable* pric);

    void display(Priceable* pric, bool quasiMC);

private:
    Priceable* priceable;  // Store the Priceable object for later use
    void displayDetails(const std::string& description, double price, double variance, const std::vector<double>& confInterval, int nameWidth, int priceWidth, int varianceWidth, int ciWidth);
};