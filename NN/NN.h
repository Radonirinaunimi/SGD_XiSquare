#pragma once

#include <string>
#include <ctime>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <eigen3/Eigen/Dense>

// Define a class Network
class Network{

	public:
		// Define the function which will define the structure of the NN
		Network(std::vector<int> nn);

		// Define the vectors which will contain the Weights and Biases
		std::vector<Eigen::MatrixXd> W;
		std::vector<Eigen::VectorXd> B;

		// Define the method that does the feed forward
		Eigen::VectorXd feed_forward(std::vector<double> input);

		// Define the method that does the Backpropagation
		void back_propagation(std::vector<double> y_exp);

	private:
		// Define the number of hidden layers
		int hidden_layers;

		// Define the vector which will contain the hidden layers
		std::vector<Eigen::VectorXd> a;

		// Define the vectors which will contain the derivatives of W & B
		std::vector<Eigen::MatrixXd> dW;
		std::vector<Eigen::VectorXd> dB;

};

