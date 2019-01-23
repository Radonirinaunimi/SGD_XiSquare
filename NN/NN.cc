#include "NN.h"


Network::NN_structure(std::vector<int> nn)
{
	// Define the number of the hidden layers
	this->hidden_layers = nn.size() - 2;

	// Initialize the vectors
	a = std::vector<Eigen::MatrixXd> (hidden_layers +  2);
	W = std::vector<Eigen::MatrixXd> (hidden_layers +  1);
	B = std::vector<Eigen::MatrixXd> (hidden_layers +  1);


	// Loop over the structure of the NN
	for (int i = 0; i < nn.size() - 1; i++)
	{
		W[i] = Eigen::MatrixXd::Random(i, i+1);
		B[i] = Eigen::MatrixXd::Random(1, i+1);
	}
}
