#include "toy.h"

using namespace std;

// Define the constructor of the class ToyPDF    
ToyPDF::ToyPDF(string const& name, double a1, double a2, int flag):

    NNPDF::PDFSet(name, 1, NNPDF::PDFSet::erType::ER_NONE)  {

    // srand(time(0));

    this -> a1 = a1;
    this -> a2 = a2; 
    this -> flag = flag;

}

// Define a Toy_NN (Basically a polynomail)
// The polynomial is of the form y(x) = a * x +b;
NNPDF::real ToyPDF::ComputePol(const vector<NNPDF::real> input) const {
    return (pow(a1, 2) * input[0] + a1 * a2);
}

// Define the method that computes the Derivative of the polynomial with respec to the 1s parameter
NNPDF::real ToyPDF::ComputePolDer1(const vector<NNPDF::real> input) const {
    return (2 * a1 * input[0] + a2);
}

// Define the method that computes the Derivative of the polynomial with respec to the 1s parameter
NNPDF::real ToyPDF::ComputePolDer2(const vector<NNPDF::real> input) const {
    return a1;
}

// Define the Numerical Derivative
NNPDF::real ToyPDF::ComputePolNum1(const vector<NNPDF::real> input) const {
	double h = 1.25e-2;
	return ( (pow(a1, 2) * input[0] + a1 * a2) - (pow((a1 - h), 2) * input[0] + (a1 - h) * a2 ) ) / h;
}

// Here, we define the method GetPDF from the class toyPDF
// Note: NNPDF::real is, basically, a float
    void ToyPDF::GetPDF(NNPDF::real const& x, NNPDF::real const& Q2, int const& n, NNPDF::real* pdf) const 
    {
        // First check the input makes sense
        int ifail = 0;
        if ( x > 1 or x < 0 ) {
            cout << "Wrong value of x passed to GetPDF: x = " << x << endl;
            ifail = 1;
        }
        if ( Q2 < 0 ) {
            cout << "Neative scale passed to GETPDF: Q2 = " << Q2 << endl;
            ifail = 1;
        }

        if (ifail) {
            cerr << "GetPDF failed with error" << endl;
            exit(-1);
        }
        

        vector<NNPDF::real> xvals{x}; 

	// Choose the appropriate GetPDF
	if (flag == 0){
		for (int i = 0; i < 14; i++){
			pdf[i] = ComputePol(xvals);
		}
	} else if (flag == 1){
		for (int i = 0; i < 14; i++){
			pdf[i] = ComputePolDer1(xvals);
		}
	} else if (flag == 2){
		for (int i = 0; i < 14; i++){
			pdf[i] = ComputePolDer2(xvals);
		}
	} else {
		for (int i = 0; i < 14; i++){
			pdf[i] = ComputePolNum1(xvals);
		}
	}

    }


// Compute the chi2 manually
NNPDF::real Chi2(int const nData, const double* data, Eigen::MatrixXd Matrix, NNPDF::real *const& theory)
{
    // Construct an array of length nData, given as input
    NNPDF::real res = 0.0;
    // Sum over the data and theory arrays
    for(int i = 0; i < nData; i++)
    {
        for(int j=0; j < nData; j++)
        {
            res += (data[i]-theory[i]) * Matrix(i,j) * (data[j]-theory[j]);
        }
    }
    return res;
}

// Compute manually the derivative of Chi2
NNPDF::real Chi2Der(int const nData, const double* data, Eigen::MatrixXd Matrix, NNPDF::real *const& theory, NNPDF::real *const& theoryDer)
{
    NNPDF::real res = 0.0;
    // Sum over the data and theory arrays
    for(int i = 0; i < nData; i++)
    {
        for(int j=0; j < nData; j++)
        {
            res += Matrix(i,j) * ( (theory[i] - data[i]) * theoryDer[j] + (theory[j] - data[j]) * theoryDer[i] );
        }
    }
    return res;
}

// Compute manually the derivative of Chi2 for DY
NNPDF::real Chi2DYDer(int const nData, const double* data, Eigen::MatrixXd Matrix, NNPDF::real *const& theory, NNPDF::real *const& theoryDer1, NNPDF::real *const& theoryDer2)
{
    NNPDF::real res = 0.0;
    // Sum over the data and theory arrays
    for(int i = 0; i < nData; i++)
    {
        for(int j=0; j < nData; j++)
        {
            res += Matrix(i,j) * ( (theory[i] - data[i]) * (theoryDer1[j] + theoryDer2[j]) + (theory[j] - data[j]) * (theoryDer1[i] + theoryDer2[i]) );
        }
    }
    return res;
}

// Invert the  Covariant Matrix
void InvertMatrix(NNPDF::matrix<double> const& inputMatrix, Eigen::MatrixXd &outputMatrix, int length)
{   
    Eigen::MatrixXd temporalMatrix(length, length);
    for(int i=0; i<length; i++)
    {
        for(int j=0; j<length; j++)
        {
            temporalMatrix(i,j) = inputMatrix(i,j);
        }
    }
    outputMatrix = temporalMatrix.inverse();
}

