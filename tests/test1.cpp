#include "../GridFunction.hpp"
#include <iostream>

const int dim = 3;

int main () {
    GridFunction<float, dim> gf1;
    gf1.set_size(GridIndex<dim>({ {10, 10, 10} }));
    gf1.minbounds = { 0, 0, 0 };
    gf1.stepsize = { 1, 1, 1 };
    for (size_t i=0; i<10; i++) 
	for (size_t j=0; j<10; j++)
	    for (size_t k=0; k<10; k++)
		gf1(GridIndex<dim>({{ i,j,k}})) = sin(i)*sin(j)*sin(k);
    std::cout << "Testing GridFunction (3D).\n";
    std::cout << "Test 1: interpolation points and weights.\n";
    GridCoordinate<dim> p = { 5.8, 5.8, 5.8 };
    auto ip = gf1.NLinearCoefficients(p);
    std::cout << "Interpolation points and weights are: \n";
    for (auto ip_it = ip.begin(); ip_it != ip.end(); ip_it++) {
	std::cout << "[";
	for (size_t i=0; i< dim; i++)
	    std::cout << ip_it->first[i] << "\t";
	std::cout << "] -> " << ip_it->second << std::endl;
    }
    std::cout << "Interpolated value is: " << gf1.NLinearValue(p) << std::endl;
    std::cout << "Precise value is: " << sin(5.8)*sin(5.8)*sin(5.8) << std::endl;
    
    return 0 ;
}

