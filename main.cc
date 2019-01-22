#include <iostream>
#include <string>
#include <stdlib.h>
#include <NNPDF/experiments.h>
#include <NNPDF/pathlib.h>
#include <NNPDF/lhapdfset.h>
#include <NNPDF/thpredictions.h>
#include <NNPDF/chisquared.h>

#include "toypdf/toy.h"


using namespace NNPDF;
using namespace std;


// This function generates the FK tables
// It takes as inputs the name of the exp. and the theory ID
FKSet LoadFK(string const& setname, int theoryid)
{
	// Initialize a vector of FK table
	vector<FKTable*> nFK;
	// Add the following element to the vector
	nFK.push_back(new FKTable{
			get_data_path() + "theory_" + to_string(theoryid) + "/fastkernel/FK_" + setname + ".dat"
			});
	// The function 'parseOperator' is located in the fkset.cc
	// It returns a value of type 'SigmaOp'??
	auto op = FKSet::parseOperator("NULL");
	return FKSet(op, nFK);
}

// Load the FK Tables
FKTable LoadFkT(string const& setname, int theoryid)
{
	FKTable fktable =  FKTable{
		get_data_path() + "theory_" + to_string(theoryid) + "/fastkernel/FK_" + setname + ".dat"
	};
	return fktable;
}

// The function takes as inputs the name of the exp. and the
// theory ID that is going to be used
DataSet LoadDataSet(string const& setname, int theoryid)
{
	// load data files and systematics treatment
	auto cd = CommonData::ReadFile(
			get_data_path() + "commondata/DATA_" + setname + ".dat",
			get_data_path() + "commondata/systypes/SYSTYPE_" + setname + "_DEFAULT.dat"
			);

	auto fk = LoadFK(setname, theoryid);

	return DataSet(cd, fk);
}



int main()
{

	string experiment;
	cout << "Type 1 for DIS and 2 for DY:" << endl;
	int input;
	cin >> input;
	while (input != 1 && input != 2){
		cout << "The value entered is not valid! Try again!" << endl;
		cin >> input;
	}
	if (input == 1){
		cout << "DIS experiment has been choosen." << endl;
		experiment = "NMC";
	} else if (input == 2){
		cout << "DY experiment has been choosen." << endl;
		experiment = "CMSJETS11";
	}

	// Allocate DIS dataset
	const auto ds = LoadDataSet(experiment, 53);
	const auto fkSet = LoadFK(experiment, 53);
	const auto fkTab = LoadFkT(experiment, 53);
	// Allocate DIS experiment (with 1 dataset)
	Experiment exp{{ds}, experiment};

	// Test DIS using the LHAPDF
	LHAPDF::setVerbosity(0);
	LHAPDFSet pdf("NNPDF31_nnlo_as_0118", PDFSet::erType::ER_NONE);
	ThPredictions LHAtheo{&pdf, &exp};
	NNPDF::real chi2lhapdf = 0;
	ComputeChi2(&exp, 1, LHAtheo.GetObs(), &chi2lhapdf);
	cout << "For DIS-LHAPDF we get: " << chi2lhapdf << endl;
	cout << "\n";

	// Using the ToyPDF
	// Initialize the coefficient of the polynomial
	// InitializePol();
	srand(time(0));
	double a1 = ((double) rand() / RAND_MAX) + 1;
	double a2 = ((double) rand() / RAND_MAX) + 1;

	cout << "***" << "\n";
	cout << "ToyPDF initialized with parameters: " << endl;
	printf("a1 = %f, a2 =%f \n", a1, a2);
	printf("The polynomial is then of the form y(x) = %f x + %f. \n", a1, a2);
	cout << "\n";

	// Initialize toypdf
	ToyPDF pdfToy("toypdf", a1, a2, 0);
	// Convolute the FK Tables with pdfs (in this case toypdf)
	ThPredictions theo{&pdf, &exp}; // Check the convolution with FKTables
	// Extract the observable from 'pdftoy'
	NNPDF::real *Obs = theo.GetObs();

	// Initialize the derivative of toypdf with respect to the 1st parameter
	ToyPDF pdfToy1("toypdf1", a1, a2, 1);
	ThPredictions theoDer1{&pdfToy1, &exp};
	// Extract the observable from 'pdftoy'
	NNPDF::real *Obs1 = theoDer1.GetObs();

	// Initialize the derivative of toypdf with respect to the 2nd parameter
	ToyPDF pdfToy2("toypdf2", a1, a2, 2);
	ThPredictions theoDer2{&pdfToy2, &exp};
	// Extract the observable from 'pdftoy'
	NNPDF::real *Obs2 = theoDer2.GetObs();

	// Initialize the derivative of toypdf with respect to the 2nd parameter
	ToyPDF pdfToyN1("toypdfN1", a1, a2, 3);
	ThPredictions theoNum1{&pdfToyN1, &exp};
	// Extract the observable from 'pdftoy'
	NNPDF::real *ObsN1 = theoNum1.GetObs();

	// Assign the experimental data into 'nData'
	const int nData = exp.GetNData();
	const double *expData = exp.GetData();
	// Get the Cholesky matrix
	const matrix<double> &covMatCholesky = exp.GetSqrtCov();
	// Get the covariance matrix
	const matrix<double> &covMat = exp.GetCovMat();
	// Invert the covariance matrix
	Eigen::MatrixXd invCovMat(nData, nData);
	InvertMatrix(covMat, invCovMat, nData);


	cout << "1. Computation of chi2: (for ToyPDF)"<< endl;

	// Compute chi_square using toypdf and the ComputeChi2_basic method in NNPDF
	NNPDF::real chi2toypdf1 = 0;
	ComputeChi2_basic(nData, 1, expData, covMatCholesky, Obs, &chi2toypdf1);

	cout << "*** <<Using the built in ComputeChi2_basic in NNPDF>>" << "\n";
	cout << "The value of chi_square using ToyPDF is: " << chi2toypdf1 << endl;
	cout << "\n";

	// Compute chi_square using Toypdf and the manual implementation of chi2
	NNPDF::real chi2toypdf2 = Chi2(nData, expData, invCovMat, Obs);

	cout << "*** <<Using the manual implementation of chi2>>" << "\n";
	cout << "The value of chi_square (should be the same as before) using ToyPDF is: " << chi2toypdf2 << endl;
	cout << "\n";


	// TODO:
	// Step 3. evaluate chi2 gradient for toy PDF.

	cout << "2. Analytical computation of the gradient descent: "<< endl;

	// Compute the derivative of chi2 with respect to the 1st parameter
	NNPDF::real chi2toypdfD1 = Chi2Der(nData, expData, invCovMat, Obs, Obs1);

	cout << "*** <<Derivative with respect to the 1st parameter>>" << "\n";
	cout << "The value of Dchi_square using ToyPDF is: " << chi2toypdfD1 << endl;
	cout << "\n";

	// Compute the derivative of chi2 with respect to the 2nd parameter
	NNPDF::real chi2toypdfD2 = Chi2Der(nData, expData, invCovMat, Obs, Obs2);

	cout << "*** <<Derivative with respect to the 2nd parameter>>" << "\n";
	cout << "The value of Dchi_square using ToyPDF is: " << chi2toypdfD2 << endl;
	cout << "\n";

	// Compute the derivative of chi2 with respect to the 1st parameter
	NNPDF::real chi2toypdfN1 = Chi2Der(nData, expData, invCovMat, Obs, ObsN1);

	cout << "*** <<Numerical derivative with respect to the 1st parameter>>" << "\n";
	cout << "The value of the 1st NDchi_square using ToyPDF is: " << chi2toypdfN1 << endl;
	cout << "\n";

	// Produce some plots
	cout << "3. Plots "<< endl;

	// Store the outputs in a txt file

	string space = "    ";
	ofstream file("../Chi2Data.txt");
	int iteration = 1000;
	double m = 1.3;
	double n = ((double) rand() / RAND_MAX) + 1;
	NNPDF::real h = 1e-4;

	ToyPDF pdfToyInit("toypdf", m, n, 0);
	ThPredictions theoInit{&pdfToyInit, &exp};
	NNPDF::real *ObsInit = theoInit.GetObs();
	NNPDF::real chi2Init = Chi2(nData, expData, invCovMat, ObsInit);

	if (input == 1){
		// Compute the DIS
		for (int i = 0; i < iteration; i++){

			// Update the value of m
			m -= (double) h;

			// Update the toypdf
			ToyPDF pdfToyNx("toypdfNx", m, n, 0);
			ThPredictions theoNumx{&pdfToyNx, &exp};
			NNPDF::real *ObsNx = theoNumx.GetObs();
			NNPDF::real chi2Nx = Chi2(nData, expData, invCovMat, ObsNx);

			// Numerical
			NNPDF::real chi2toypdfNx = (chi2Init - chi2Nx)/h;
			chi2Init = chi2Nx;

			// Analytical
			ToyPDF pdfToyDx("toypdfDx", m, n, 1);
			ThPredictions theoDerx{&pdfToyDx, &exp};
			NNPDF::real *ObsDx = theoDerx.GetObs();
			NNPDF::real chi2toypdfDx = Chi2Der(nData, expData, invCovMat, ObsNx, ObsDx);

			// file << n << "  " << chi2toypdfDx << "  " << chi2toypdfNx << endl;
			file.width(10);
			file << m << "  ";
			file.width(10);
			file << chi2toypdfDx << "  ";
			file.width(10);
			file << chi2toypdfNx << "  ";
			file.width(10);
			file << chi2toypdfNx/chi2toypdfDx << "\n";
		}
		file.close();
	} else if (input == 2){
		// Compute the DY
		for (int i = 0; i < iteration; i++){

			// Update the value of m
			m -= (double) h;

			// Update the toypdf
			ToyPDF pdfToyNx("toypdfNx", m, n, 0);
			ThPredictions theoNumx{&pdfToyNx, &exp};
			NNPDF::real *ObsNx = theoNumx.GetObs();
			NNPDF::real chi2Nx = Chi2(nData, expData, invCovMat, ObsNx);

			// Numerical
			NNPDF::real chi2toypdfNx = (chi2Init - chi2Nx)/h;
			chi2Init = chi2Nx;

			// Analytical
			ToyPDF pdfToyDx("toypdfDx", m, n, 1);
			// First
			ThPredictions theoDerx1{&pdfToyDx, &pdfToyNx, &fkTab};
			NNPDF::real *ObsDx1 = theoDerx1.GetObs();
			// Second
			ThPredictions theoDerx2{&pdfToyNx, &pdfToyDx, &fkTab};
			NNPDF::real *ObsDx2 = theoDerx2.GetObs();
			// Combined
			NNPDF::real chi2toypdfDx = Chi2DYDer(nData, expData, invCovMat, ObsNx, ObsDx1, ObsDx2);

			// file << n << "  " << chi2toypdfDx << "  " << chi2toypdfNx << endl;
			file.width(10);
			file << m << "  ";
			file.width(10);
			file << chi2toypdfDx << "  ";
			file.width(10);
			file << chi2toypdfNx << "  ";
			file.width(10);
			file << chi2toypdfNx/chi2toypdfDx << "\n";
		}
		file.close();
	}


	return 0;
}
