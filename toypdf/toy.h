#pragma once

#include <string>
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <eigen3/Eigen/Dense>

#include <NNPDF/experiments.h>
#include <NNPDF/pathlib.h>
#include <NNPDF/lhapdfset.h>
#include <NNPDF/thpredictions.h>
#include <NNPDF/pdfset.h>

// Define a class ToyPDF which inherites from the
// main class NNPDF
class ToyPDF: public NNPDF::PDFSet {

    private:
	// Define the method that caomputes the NN
        NNPDF::real ComputePol(const std::vector<NNPDF::real>) const;
        NNPDF::real ComputePolDer1(const std::vector<NNPDF::real>) const;
        NNPDF::real ComputePolDer2(const std::vector<NNPDF::real>) const;
        NNPDF::real ComputePolNum1(const std::vector<NNPDF::real>) const;
        // Define the method that initializes the generation of the random numbers
        double a1, a2;
	int flag;

    public:
	// Define the constructor that is going to be called
	// each time we call ToyPDF
        ToyPDF(std::string const&, double, double, int);
	// Define the method that gets the PDF
        void GetPDF(NNPDF::real const& x, NNPDF::real const& Q2, int const& n, NNPDF::real* pdf) const;

};

// Define the computation of chi_square
NNPDF::real Chi2(int const nData, const double* data, Eigen::MatrixXd Matrix, NNPDF::real *const& theory);
// Define the computation of chi2 for the derivative with respect to the 1st/2nd parameter
NNPDF::real Chi2Der(int const nData, const double* data, Eigen::MatrixXd Matrix, NNPDF::real *const& theory, NNPDF::real *const& theoryDer);
// Define the computation of chi2 for the derivative with respect to the 1st/2nd parameter
NNPDF::real Chi2DYDer(int const nData, const double* data, Eigen::MatrixXd Matrix, NNPDF::real *const& theory, NNPDF::real *const& theoryDer1, NNPDF::real *const& theoryDer2);
// Define a function which inverts the Matrix
void InvertMatrix(NNPDF::matrix<double> const&, Eigen::MatrixXd &, int length);

