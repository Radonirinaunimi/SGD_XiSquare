#include<iostream>
#include<vector>
#include<cmath>
#include<eigen3/Eigen/Dense>
#include<eigen3/Eigen/Core>

// Define the Sigmoid function
double Sigmoid(double x)
{
	return 1 / (1 + std::exp(-x));
}

std::vector<double> _dot0(Eigen::MatrixXd M, Eigen::MatrixXd v){

	std::vector<double> result;
	double value;

	for(int i = 0; i < M.rows(); i++){
		value = (M.row(i).array() * v.array()).sum();
		result.push_back(value);
	}

	return result;

}

Eigen::VectorXd _dot(Eigen::MatrixXd M, Eigen::VectorXd v){

	std::vector<double> result;
	double value;

	for(int i = 0; i < M.rows(); i++){
		value = (M.row(i).array() * v.transpose().array()).sum();
		result.push_back(value);
	}

	Eigen::VectorXd test = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned> (result.data(), result.size());
	return test;

}

int main(){

	// Define a Matrix with undefined size
	Eigen::MatrixXd M1(3, 3);
	// Fill the Matrix with some numbers
	M1 << 1, 2, 3,
	      4, 5, 6,
	      7, 8, 9;
	// std::cout << M1 << "\n";
	Eigen::MatrixXd Vm = Eigen::MatrixXd::Random(1, 3);
	std::cout << "Vector in Matrix: " << Vm << "\n";

	// Print each row of M1
	std::cout << "\n" << std::endl;
	for(int i = 0; i < M1.rows(); i++){
		std::cout << "Row " << i << " of M1: " << M1.row(i) << "\n";
		std::cout << "Row " << i << " of M1.array() square: " << M1.row(i).array() * M1.row(i).array() << "\n";
		std::cout << "Sum of the element: " << (M1.row(i).array() * M1.row(i).array()).sum() << "\n";
		// Test with vector
		std::cout << " Transform Vec to array: " << M1.row(i).array() * Vm.array() << "\n";
	}

	// Define a matrix filled with random numbers of type double
	std::cout << "\n" << std::endl;
	Eigen::MatrixXd M2 = Eigen::MatrixXd::Random(3, 3);
	// std::cout << M2 << "\n";

	// Define a vector of Matrices
	std::vector<Eigen::MatrixXd> v_M;

	for(int i = 0; i < 4; i++){
		v_M.push_back(Eigen::MatrixXd::Random(3, 3));
	}

	for(int i = 0; i <v_M.size(); i++){
		std::cout << v_M[i] << std::endl;
		std::cout << "\n" << std::endl;
	}

	// Do an element-wise product between matrices
	std::cout << "\n" << std::endl;
	std::cout << "Elementwise product between Matrices: \n" << M1.array() * M1.array() << "\n";

	// Test the function which computes the _dot product
	std::vector<double> vX1 = _dot0(M1, Vm);
	std::cout << "\n" << std::endl;
	std::cout << "Result: " << "\n";
	for(int i = 0; i <vX1.size(); i++){
		std::cout << vX1[i] << std::endl;
	}


	Eigen::VectorXd Vm2 = Eigen::VectorXd::Random(3);
	std::cout << "This is the result: " << (M1.row(1).array() * Vm2.transpose().array()).sum()  << "\n";
	Eigen::VectorXd vX2 = _dot(M1, Vm2);
	std::cout << "\n" << std::endl;
	std::cout << "Comparison to the above result: \n" << vX2 << std::endl;

	// Apply the Sigmoid to Eigen elements
	std::cout << "\n" << std::endl;
	std::cout << "When applying the Sigmoid: " << "\n";
	std::cout << vX2.unaryExpr(&Sigmoid) << std::endl;


	return 0;

}
