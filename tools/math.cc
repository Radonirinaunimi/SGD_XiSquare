#include "math.h"

// Define the Sigmoid function
double Sigmoid(double x)
{
        return 1 / (1 + std::exp(-x));
}


// Define the dot product
Eigen::VectorXd _dot(Eigen::MatrixXd M, Eigen::VectorXd v)
{

        std::vector<double> result;
        double value;

        for(int i = 0; i < M.rows(); i++)
	{
                value = (M.row(i).array() * v.transpose().array()).sum();
                result.push_back(value);
        }

        Eigen::VectorXd vec = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned> (result.data(), result        .size());
        return vec;
}


// Define the Hadamard product
Eigen::VectorXd Hadamard(Eigen::VectorXd v1, Eigen::VectorXd v2)
{
        std::vector<double> result;
        double value;

        for(int i = 0; i < v1.size(); i++)
        {
                value = v1(i) * v2(i);
                result.push_back(value);
        }


        Eigen::VectorXd test = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned> (result.data(), result.size());
        return test;
}

