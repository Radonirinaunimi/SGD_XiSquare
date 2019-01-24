#include "NN.h"


Network::Network(std::vector<int> nn)
{
	// Define the number of the hidden layers
	this->hidden_layers = nn.size() - 2;

	// Initialize the vectors
	a = std::vector<Eigen::VectorXd> (hidden_layers +  2);
	W = std::vector<Eigen::MatrixXd> (hidden_layers +  1);
	B = std::vector<Eigen::VectorXd> (hidden_layers +  1);


	// Loop over the structure of the NN
	for (int i = 0; i < hidden_layers + 1; i++)
	{
		W[i] = Eigen::MatrixXd::Random(nn[i], nn[i+1]);
		B[i] = Eigen::VectorXd::Random(nn[i+1]);
	}

	std::cout << "It works!!" << std::endl;

	// Try to print some results
        for(int i = 0; i < B.size(); i++){
                std::cout << B[i] << std::endl;
                std::cout << "\n" << std::endl;
        }

}


// Eigen::MatrixXd Neural::feed_forward(std::vector<double> input)
// {
// 	// Define the first layer to be the input
// 	a[0] = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned> (input.data(), input.size());
//
// 	// Compute the next layer
// 	for (int i = 1; i < hidden_layers + 2; i++)
// 	{
// 		a[i] = a[i-1];
// 	}
// }
