#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
 	/**
  	* TODO:
    * Calculate the RMSE here.
  	*/
	VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
        // Ignore data if  it's not complete
        cout << " CalculateRMSE : Invalid estimations or  ground_truth vector size" << endl;
        return rmse;
    }
    
    for (int idx = 0; idx < estimations.size(); ++idx) {
        Eigen::VectorXd diff = estimations[idx] - ground_truth[idx];
        diff  = diff.array() * diff.array();
        rmse += diff;
    }
    
    //calculate the mean
    rmse = rmse / estimations.size();
    
    //calculate the squared root
    rmse = rmse.array().sqrt();
    
    return rmse;
}