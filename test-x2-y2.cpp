#include <iostream>
#include <vector>
#include <memory>
#include "nelder-mead.h"

// f(x,y) = (x-6.6)^2 + (y-18)^2
double paraboloid(std::vector<double>& input)
{
  double x = input[0];
  double y = input[1];
  
  return (x-6.6)*(x-6.6) + (y-18)*(y-18);
}

int main()
{
 std::vector<double> start{-12.0, 65.0};
 std::vector<double> output;
 NelderMead<double> minimizer(start, paraboloid);
 minimizer.setError(1e-8);
 minimizer.minimize(output);
 minimizer.printSimplex();

 std::cout << std::endl << "paraboloid output minimum vector [x,y]:" << std::endl;
 for (auto& vertex: output)
   std::cout << vertex << " ";
 std::cout << std::endl;

return 0;
}
