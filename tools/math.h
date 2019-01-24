#include<iostream>
#include<vector>
#include<cmath>
#include<eigen3/Eigen/Dense>
#include<eigen3/Eigen/Core>


// Initialize the Sigmoid function
double Sigmoid(double x);

// Initialize the dot product
Eigen::VectorXd _dot(Eigen::MatrixXd M, Eigen::VectorXd v);
