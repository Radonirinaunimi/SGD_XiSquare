#pragma once

#include <string>
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <eigen3/Eigen/Dense>

// Define a class Network
class Network{

	public:
		// Define the function which will define the structure of the NN
		NN_structure(std::vector<int> nn);

		// Define the vectors which will contain the Weights and Biases
		std::vector<Eigen::MatrixXd> W;
		std::vector<Eigen::MatrixXd> B;

		// Define the method that does the feed forward
		Eigen::MatrixXd feed_forward(std::vector<double> input);

		// Define the method that does the Backpropagation
		void back_propagation(std::vector<double> y_exp);

	private:
		// Define the number of hidden layers
		int hidden_layers;

		// Define the vector which will contain the hidden layers
		std::vector<Eigen::MatrixXd> a;

		// Define the vectors which will contain the derivatives of W & B
		std::vector<Eigen::MatrixXd> dW;
		std::vector<Eigen::MatrixXd> dB;

};
