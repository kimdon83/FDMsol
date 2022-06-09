#ifndef _VECTOR_H_SHJ_
#define _VECTOR_H_SHJ_

#include "define.h"
#include <math.h>

class CSparceMatrix;

class CVector {
public:
	int dim; //dimension
	real *value; //elements

public:
	CVector(void) //Constructor
	{
		dim = 0;
		value = NULL;
	}
	CVector(int d) //Constructor with a dimension
	{
		SetDimension(d);
	}
	~CVector(void) //Destructor
	{
		delete [] value;
	}

public:
	//Set the dimension
	void SetDimension(int d)
	{
		int i;

		dim = d;
		value = new real [dim];
		for(i=0; i<dim; i++) value[i] = 1./(real)dim;
	}
	int GetDimension(void) { return dim; }

	//Set a value
	void SetValue(int loc, real val) { value[loc] = val; }
	//Set values
	void SetValues(real *val) { memcpy(value, val, sizeof(real)*dim); }
    void SetValues(CVector &Vec) { memcpy(value, Vec.value, sizeof(real)*dim); }

	//Get a value
	real GetValue(int loc) { return value[loc]; }
    real* GetValues(void)  { return value; }

    //Return the maximum relative difference
    real GetMaxRelDiff(CVector &Vec)
    {
        int i;
        real maxErr = 0., err;
#pragma omp parallel for private(err) shared(maxErr)
        for(i=0; i<dim; i++) {
            err = fabs((value[i]-Vec.value[i])/Vec.value[i]);
            if(err > maxErr) maxErr = err;
        }

        return maxErr;
    }

public:
    //this = this - mul*Vec
    real SubtractAMultipliedVector(real mul, CVector &Vec)
    {
        int i;
        real sum=0.;
        for(i=0; i<dim; i++) {
            value[i] -= mul*Vec.value[i];
            sum += value[i];
        }
        return sum;
    }

	//Perform an inner-production of a matrix and a vector and store it
	real Multiply(CSparceMatrix &Mat, CVector &Vec);
    //Perform a scalar multiplication
    void Multiply(real c)
    {
        int i;
        for(i=0; i<dim; i++) value[i] *= c;
    }
    //Perform a Vec^T.this
    real Multiply(CVector &Vec)
    {
		int i;
		real sum = 0.;
#pragma omp parallel for private(i) reduction(+:sum)
		for (i = 0; i < dim; i++) {
			sum += this->value[i] * Vec.value[i];
		}
		return sum;
    }
	//Solve the equation of Mat.this = Vec
	void SolveByPointJacobiIterativeMethod(CSparceMatrix &Mat, CVector &Vec);
    void SolveByGaussSeidelIterativeMethod(CSparceMatrix &Mat, CVector &Vec);
    void SolveBySuccessiveOverRelaxationMethod(CSparceMatrix &Mat, CVector &Vec, real omega);

	//Print the vector
	void Print(ostream &outs)
	{
		int i;

        outs.precision(5);
		outs.setf(ios::fixed, ios::floatfield);
		for(i=0; i<dim; i++) {
			outs << setw(14) << value[i] << ",";// << endl;

		}
	}
};

#endif