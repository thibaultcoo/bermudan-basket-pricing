#pragma once

#include "instantiating.h"
#include "pricing.h"

// Display object for the pricing results
class Display
{
    public:
        // Constructor to initialize with a Priceable object
        explicit Display(Priceable* pric);

        // Function that prompts out our pricing results in a windows terminal
        void display(Priceable* pric, bool quasiMC);

    private:
        // Helper function to prettify our output table
        void displayDetails(const std::string& description, double price, double variance, const std::vector<double>& confInterval, int nameWidth, int priceWidth, int varianceWidth, int ciWidth);
};

// Display object for the initial parameters
class paramsDisplay {
    public:
        // Constructor to initialize display object with a const reference to an Instance
        explicit paramsDisplay(const Instance& inst) {}

        // Function that will prompt out the parameters used for this instance
        void display(const Instance& inst);
};