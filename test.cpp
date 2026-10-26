#include <iostream>
#include <vector>
#include <memory>
#include "nelder-mead.h"

// rosenbrock banana function
double banana(const std::vector<double>& input)
{
  // x = input[0]
  // y = input[1]
  // f(x,y) =  (1-x)^2 + 100(y-x^2)^2	
  return (1.0-input[0])*(1.0-input[0]) + 100*(input[1]-input[0]*input[0])*(input[1]-input[0]*input[0]);
}

int main()
{
 std::vector<double> start{12.0, -5.0};
 std::vector<double> output;
 NelderMead<double> minimizer(start, banana);
 minimizer.minimize(output);
 minimizer.printSimplex();

 std::cout << std::endl << "output minimum vector [x,y]:" << std::endl;
 for (auto& vertex: output)
   std::cout << vertex << " ";
 std::cout << std::endl;

return 0;
}
