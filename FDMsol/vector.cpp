#include "define.h"
#include "vector.h"
#include "sparce_matrix.h"


#include <math.h>

real CVector::Multiply(CSparceMatrix &Mat, CVector &Vec)
{
	int i, j;
	real sum = 0.;

	//Check dimensions
	if ((dim != Mat.GetDimension()) || (dim != Vec.GetDimension())) {
		cout << "Dimension mismatch at CVector::Multiply" << endl;
		exit(0);
	}

	//Initialize the elements
	memset(value, 0, sizeof(real)*dim);
#pragma omp parallel for private(j) 
	//Perform the inner-production
	for (i = 0; i < dim; i++)
	{
		for (j = 0; j < Mat.nStripe; j++) {
			if (Mat.locElem[i][j] >= 0)
				value[i] += Mat.elem[i][j] * Vec.value[Mat.locElem[i][j]];
		}
	}
	sum = 0;
	for (i = 0; i < dim; i++) sum += value[i];

	return sum;
}

//Solve the equation of Mat.this = Vec by the Point Jacobi Iterative Method
void CVector::SolveByPointJacobiIterativeMethod(CSparceMatrix &Mat, CVector &Vec)
{
	int i, j;
	real term, mul;
	CVector PrevFlux, Err;
	real maxErr = 0., relErr;
	int nIt = 0;

	//Check dimensions
	if ((dim != Mat.GetDimension()) || (dim != Vec.GetDimension())) {
		cout << "Dimension mismatch at CVector::Multiply" << endl;
		exit(0);
	}

	//Allocate the previous flux vector
	PrevFlux.SetDimension(dim);
	Err.SetDimension(dim);


	do {
		PrevFlux.SetValues(this->value);
		//Update flux
		maxErr = 0.;

		#pragma omp parallel for private(j,term,mul)
		for (i = 0; i < dim; i++)
		{
			term = Vec.value[i]; mul = 0.;
			for (j = 0; j < Mat.nStripe; j++) {
				if (Mat.locElem[i][j] >= 0) {
					if (Mat.locElem[i][j] == i) mul = Mat.elem[i][j];
					else term -= Mat.elem[i][j] * PrevFlux.value[Mat.locElem[i][j]];
				}
			}
			if (mul == 0.) {
				cout << "A diagonal term of the matrix is zero at CVector::SolveByPointJacobiIterativeMethod" << endl;
				exit(0);
			}
			this->value[i] = term / mul;

			Err.value[i] = this->value[i] - PrevFlux.value[i];
		}

		nIt++;

		relErr = Err.Multiply(Err) / this->Multiply(*this);
		//cout << relErr << endl;
	} while (relErr > 1.e-6);
	cout << "Inner Iter=" << nIt;
}

//Solve the equation of Mat.this = Vec by the Gauss Sedel Iterative Method
void CVector::SolveByGaussSeidelIterativeMethod(CSparceMatrix &Mat, CVector &Vec)
{
	int i, j;
	real term, mul;
	CVector PrevFlux, Err;
	real maxErr = 0., relErr;
	int nIt = 0;


	//Check dimensions
	if ((dim != Mat.GetDimension()) || (dim != Vec.GetDimension())) {
		cout << "Dimension mismatch at CVector::Multiply" << endl;
		exit(0);
	}

	//Allocate the previous flux vector
	PrevFlux.SetDimension(dim);
	Err.SetDimension(dim);

	do {
		//Solve the matrix equation
		PrevFlux.SetValues(this->value);
		//Update flux
		maxErr = 0.;

#pragma omp parallel for private(j,term,mul)

		for (i = 0; i < dim; i++)
		{
			term = Vec.value[i]; mul = 0.;
			for (j = 0; j < Mat.nStripe; j++) {
				if (Mat.locElem[i][j] >= 0) {
					if (Mat.locElem[i][j] == i) mul = Mat.elem[i][j];
					else term -= Mat.elem[i][j] * this->value[Mat.locElem[i][j]];
				}
			}
			if (mul == 0.) {
				cout << "A diagonal term of the matrix is zero at CVector::SolveByGaussSeidelIterativeMethod" << endl;
				exit(0);
			}
			this->value[i] = term / mul;

			Err.value[i] = this->value[i] - PrevFlux.value[i];
		}
		nIt++;

		relErr = Err.Multiply(Err) / this->Multiply(*this);
	} while (relErr > 1.e-6);
	cout << "Inner Iter=" << nIt;

}

//Solve the equation of Mat.this = Vec by the Successive-Over-Relaxation Method
void CVector::SolveBySuccessiveOverRelaxationMethod(CSparceMatrix &Mat, CVector &Vec, real omega)
{
	int nIt = 0;

	int i, j;
	real term, mul;
	CVector PrevFlux, Err;
	real maxErr = 0., relErr;

	//Check dimensions
	if ((dim != Mat.GetDimension()) || (dim != Vec.GetDimension())) {
		cout << "Dimension mismatch at CVector::Multiply" << endl;
		exit(0);
	}

	//Allocate the previous flux vector
	PrevFlux.SetDimension(dim);
	Err.SetDimension(dim);

	do {
		//Solve the matrix equation
		PrevFlux.SetValues(this->value);
		//Update flux
		maxErr = 0.;


#pragma omp parallel for private(j,term,mul)

		for (i = 0; i < dim; i++)
		{

			term = Vec.value[i]; mul = 0.;
			for (j = 0; j < Mat.nStripe; j++) {
				if (Mat.locElem[i][j] >= 0) {
					if (Mat.locElem[i][j] == i) mul = Mat.elem[i][j];
					else term -= Mat.elem[i][j] * this->value[Mat.locElem[i][j]];
				}
			}
			if (mul == 0.) {
				cout << "A diagonal term of the matrix is zero at CVector::SolveBySuccessiveOverRelaxationMethod" << endl;
				exit(0);
			}
			this->value[i] = term / mul;
			this->value[i] = omega*this->value[i] + (1 - omega)*PrevFlux.value[i];

			Err.value[i] = this->value[i] - PrevFlux.value[i];
		}

		nIt++;

		relErr = Err.Multiply(Err) / this->Multiply(*this);
	} while (relErr > 1.e-6);
	cout << "Inner Iter= " << nIt;

}