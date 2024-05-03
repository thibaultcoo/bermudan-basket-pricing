#include "usefulfunctions.h"

//Print vectors
void print_vect(const std::vector<double>& v)
{
    std::ostream_iterator<double> out_it(std::cout, ", ");
    std::copy(v.begin(), v.end(), out_it);
    std::cout << std::endl;
}

//Vectors to csv
void exportVectortoXl(std::string file, std::vector<double> V)
{
    std::ofstream myfile;
    myfile.open(file, std::ofstream::app);

    myLong vsize = V.size();
    for (myLong n = 0; n < vsize; n++)
    {
        myfile << V[n] << "\n";

    }
    myfile.close();
}