#include <functional>
#include <iostream>
#include <concepts>
#include <type_traits>
#include <vector>
#include <cmath>


template<std::floating_point T>
using MinFunc = std::function<T(const std::vector<T>&)>;

// Nelder-Mead class
template<std::floating_point T>
class NelderMead
{
     public:	
         NelderMead(const std::vector<T>& initial, MinFunc<T> minimize_function) : initial(initial), error(1e-6), vertex_length((T)1.0), simplex(initial.size()+1, std::vector<T>(initial.size(), 0.0)), minifunc(minimize_function) {};

	 void setVertexLen(T vertex_length);
	 void setError(double error);
         void minimize(std::vector<T>& minimum);
	 double minimizeStep();
	 void printSimplex();


     private:

         void fillSimplex();

         double error;
	 T vertex_length;
	 std::vector<T> initial;
         std::vector<std::vector<T>> simplex;
         MinFunc<T> minifunc;
};


// ---------------------------------------------------- implementation ------------------------------------------------------------


// change simplex edge size
template<std::floating_point T>
void NelderMead<T>::setVertexLen(T vertex_length) 
{
  this->vertex_length=vertex_length;
}


// set calculation error
template<std::floating_point T>
void NelderMead<T>::setError(double error) 
{
  this->error=error;
}



// create N+1 simplex points, for 2d create triangle, for 3d create tetraeder
template<std::floating_point T>
void NelderMead<T>::fillSimplex()
{
   for (size_t i=0;i<initial.size();i++)
     simplex[0][i] = initial[i];

   for (size_t i=0;i<initial.size();i++)
   {
     for (size_t j=0;j<initial.size();j++)
       if (i==j)
         simplex[i+1][j] = initial[j]+vertex_length;
       else
         simplex[i+1][j] = initial[j];

   }
}



// print simplex vertices
template<std::floating_point T>
void NelderMead<T>::printSimplex()
{
   for (auto& vertex: simplex)
   {
     std::cout << " point: ";

     for (auto& point: vertex)
       std::cout << point << " ";
     std::cout << "value: " << minifunc(vertex) << std::endl;
   }
}


// search the minimum of the function so the simplex size error is under the limit
template<std::floating_point T>
void NelderMead<T>::minimize(std::vector<T>& minimum)
{
  fillSimplex();
  
  int ssize = simplex.size();
  
  #ifdef DEBUG
    std::cout << "minimize start, simplex size:" << ssize << std::endl;
  #endif  

  double calc_error = 100.0;

  int loop_count=0;

  while(calc_error > error && loop_count<1000)
  {
    calc_error = minimizeStep();
    loop_count++;
  
    #ifdef DEBUG
      std::cout << "step: " << loop_count << std::endl;
      std::cout << "======================================================" << std::endl;
    #endif  
  }


  std::vector<T> last = simplex[ssize-1];

  minimum.clear(); 
  for (auto& point: last) 
  {
    #ifdef DEBUG
      std::cout << point << " ";
    #endif  

    minimum.push_back(point);	
  }
}


// one step of the minimization
template<std::floating_point T>
double NelderMead<T>::minimizeStep()
{
  // largest value vertex
  size_t largest_index = 0;
  T largest_value = minifunc(simplex[0]);
   
  #ifdef DEBUG
    std::cout << "__________________________________________" << std::endl;
    std::cout << "minimize step: " << std::endl;
    std::cout << "[0] func value: " << largest_value << std::endl;
  #endif

  for(size_t i=1;i<simplex.size();i++)
  {
    T value = minifunc(simplex[i]);

    #ifdef DEBUG
      std::cout << "[" << i << "] func value: " << value << std::endl;
    #endif

    if (value > largest_value)
    {	    
      largest_index = i;
      largest_value = value;
    }      
  }
    
  #ifdef DEBUG
      std::cout << "largest index: " << largest_index << " value:" << largest_value << std::endl;
  #endif

  // reflection point
  std::vector<T> reflection(simplex[0].size(), 0);

  for (size_t i=0;i<simplex.size();i++)
  {
     if (i==largest_index)
	     continue;


     for(size_t j=0;j<reflection.size();j++)
     {
       reflection[j] += simplex[i][j];

       #ifdef DEBUG
         std::cout << simplex[i][j] << " ";
       #endif
     }
     
     #ifdef DEBUG
       std::cout << std::endl;
     #endif
  }	  
  
  for(auto& v : reflection)
    v /= (T)simplex[0].size();

  #ifdef DEBUG
    std::cout << "reflection point output: " << std::endl;
    for(auto& v : reflection)
      std::cout << v << " " << std::endl;
  #endif

  // reflect largest point in simplex over reflection point
  std::vector<T> newPoint(simplex[0].size(), 0);
  std::vector<T> newPointLonger(simplex[0].size(), 0);
  std::vector<T> newPointContracted(simplex[0].size(), 0);
  T newValue;
  T newValueLonger;

  for (size_t i=0;i<simplex[largest_index].size();i++)
    newPoint[i]=2*reflection[i]-simplex[largest_index][i];

  newValue = minifunc(newPoint);
  
  #ifdef DEBUG
    std::cout << "new point after reflection:" << std::endl;
    for (size_t i=0;i<newPoint.size();i++)
      std::cout << newPoint[i] << " ";
    std::cout << " value: " << newValue << std::endl;    
  #endif

  // if new value is smaller, prolong reflection, if the value is larger, then contract
  newValueLonger = 10 * newValue;

  if (newValue < largest_value)
  {	  
    for (size_t i=0;i<simplex[largest_index].size();i++)
      newPointLonger[i]=3*reflection[i]-2*simplex[largest_index][i];

    newValueLonger = minifunc(newPointLonger);
  
    #ifdef DEBUG
      std::cout << "new point after longer reflection:" << std::endl;
      for (size_t i=0;i<newPointLonger.size();i++)
        std::cout << newPointLonger[i] << " ";
      std::cout << " value: " << newValueLonger << std::endl;    
    #endif
  }

  // if new longer value is smaller use the point
  if (newValue >= newValueLonger && newValueLonger < largest_value)
  {
    for (size_t i=0;i<simplex[largest_index].size();i++)
      simplex[largest_index][i]=newPointLonger[i];
    
    #ifdef DEBUG
      std::cout << "selected new longer point, value:" << newValueLonger << std::endl;	  
    #endif
  } 
  else if (newValue < newValueLonger && newValue < largest_value)
  {
    for (size_t i=0;i<simplex[largest_index].size();i++)
      simplex[largest_index][i]=newPoint[i];
   
    #ifdef DEBUG
      std::cout << "selected new point, value: " << newValue << std::endl;	  
    #endif
  }
  else if (newValue > largest_value && newValueLonger > largest_value)
  {
    // contraction in case new value and new value longer are not smaller
    for (size_t i=0;i<simplex[largest_index].size();i++)
      newPointContracted[i]=(reflection[i]+simplex[largest_index][i])/2;

    for (size_t i=0;i<simplex[largest_index].size();i++)
      simplex[largest_index][i]=newPointContracted[i];
    
    #ifdef DEBUG
      std::cout << "selected new contracted point, value:" << newValueContracted << std::endl;	  
    #endif
  }

  #ifdef DEBUG
  printSimplex();
  #endif

  // calculate error
  #ifdef DEBUG
    std::cout << "calc error start " << std::endl;	  
  #endif

  double calc_error = 0.0;
  for (size_t i=0;i<simplex.size()-1;i++)
    calc_error += fabs(minifunc(simplex[i]) - minifunc(simplex[i+1]));

  #ifdef DEBUG
    std::cout << "calc error1: " << calc_error << std::endl;	  
  #endif	 
  calc_error /= simplex[0].size()*(1+simplex[0].size()); 


  #ifdef DEBUG
    std::cout << "calc error: " << calc_error << std::endl;	  
  #endif	 

  return calc_error; 
}


