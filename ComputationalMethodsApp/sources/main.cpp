#include "Functions.h"
#include "Vector.h"
#include "Matrix.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>

using namespace std;

int main ()
{
	ifstream F;
	F.open("../Data/input.txt");
	if (!F.is_open())
	{
		cout << "File is not exist" << endl;
		return -1;
	}

	Matrix A;
	F>>A;
	Vector b;
	//F>>b;
	F.close();
	Matrix q, r;
	double eps = 0.000001;
	//QRDecomposition(A, q, r);
	//Matrix B(A);
	//A.Identity();
	int k;

	//cout<<A<<endl;
	//cout<<A.Inverse()<<endl;

	//cout<<A.Segment()<<endl;
	//cout<<A.MaxValue(eps)<<endl;
	//cout<<Functions::GetValues(A,eps);
	cout<<Functions::BisectionMethod(A, eps);


	/*cout<<sqrt((A * A).MaxValue(eps));
	cout<<Functions::MinResidual(A, b, eps, k)<<endl;
	cout<<"MinResidual. Number of iterations: "<<endl<<k<<endl;
	
	cout<<endl<<A<<endl;

	cout<<Functions::SimpleIteration(A, b, eps, k)<<endl;
	cout<<"SimpleIteration. Theoretical number of iterations: "<<k<<endl;
	double x = ((A.MaxValue(eps) - A.MinValue(eps)) / (A.MaxValue(eps) + A.MinValue(eps)));

	cout<<"SimpleIteration. Number of iterations: "<<-log(2 / eps) / -log(1 / x)<<endl;
	cout<<endl<<A<<endl;*/

	//cout<<log(1 / x)<<endl;
	//double detA = 0.0;
	//cout<<Functions::GaussS(A, b, detA)<<endl;
	//cout<<endl<<detA<<endl<<endl;

	//cout<<A.Norm()<<endl<<b<<endl;
	//cout<<A.MinValue(0.001)<<endl;
	//cout<<A.ConditionNumber()<<endl;
	//cout<<A.Norm()<<endl;
	//cout<<A.Inverse()<<endl;
	//cout<<endl<<A;
	//cout<<A.ConditionNumber()<<endl;

	//cout<<A.Inverse();
	//cout<<A * A.Inverse();

	//cout<<Functions::RotationMethod(A, b)<<endl;


	/*cout<<detA<<endl;
	cout<<Functions::GaussS(B, b, detA)<<endl;
	
	//A.SwapLine(0,2);
	cout<<A<<endl<<B<<endl;
	cout<<A.Inverse();*/

	return 0;
}