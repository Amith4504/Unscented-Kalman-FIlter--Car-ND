#include "tools.h"
using Eigen::VectorXd;
using std::vector;
#include <iostream>

Tools::Tools() {}

Tools::~Tools() {}


VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,const vector<VectorXd> &ground_truth)

{
  /**
   * TODO: Calculate the RMSE here.
   */
  // Assertion

  VectorXd rmse = VectorXd(4);
  rmse << 0,0,0,0;
  int size  = estimations.size();


  if(size == 0 || size <= 3)
  {
    std::cout << "Estimation vector size out of bounds or less than required" << std::endl;
    return rmse;
  }
  else if(estimations.size() != ground_truth.size())
  {
    std::cout << "Dimensions do not match !!" << std::endl;
    return rmse;
  }
  else
  {
     int  n = estimations.size();

     for(int i = 0 ; i< n ; ++i)
     {
       VectorXd diff = estimations[i] - ground_truth[i];
       diff = diff.array() * diff.array();
       rmse += diff ;
     }

     rmse = rmse/n;
     rmse = rmse.array().sqrt();
     return rmse;



  }

}


