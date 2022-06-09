#ifndef _SPARCE_MATRIX_H_SHJ_
#define _SPARCE_MATRIX_H_SHJ_

#include "define.h"

class CVector;

class CSparceMatrix { 
//Declare friend classes
friend class CVector;

private:
	int dim; // dimension of the square sparce matrix
	int nStripe; //maximum number of stripes
	int **locElem; //colume indices
	               //[i][j], i = 0, ..., dim-1 (row index)
	               //        j = 0, ..., nStripe-1
	real **elem; //element
	             //[i][j], i = 0, ..., dim-1 (row index)
	             //        j = 0, ..., nStripe-1

public:
	CSparceMatrix(void) //Constructor
	{
		dim = 0; nStripe = 0;
		locElem = NULL; elem = NULL;
	}
	CSparceMatrix(int d, int s) //Constructor with a dimension and a number of stripes
	{
		SetDimension(d, s);
	}
	~CSparceMatrix(void) //Destructor
	{
		int i;

		if(locElem) {
			for(i=0; i<dim; i++) delete [] locElem[i];
			delete [] locElem;
		}
		if(elem) {
			for(i=0; i<dim; i++) delete [] elem[i];
			delete [] elem;
		}
	}

	int GetDimension(void) { return dim; }
	int GetNumberOfStripes(void) { return nStripe; }

public:
	//Set the dimension and the number of stipes
	void SetDimension(int d, int s)
	{
		int i,j;

		dim = d; nStripe = s;
		locElem = new int * [dim];
		elem = new real * [dim];
		for(i=0; i<dim; i++)
		{
			locElem[i] = new int [nStripe]; 
			for(j=0; j<nStripe; j++) locElem[i][j] = -1;
			elem[i] = new real [nStripe];
			memset(elem[i], 0, sizeof(real)*nStripe);
		}
	}

	//Set an element
	int SetElement(int row, int col, real val)
	{
		int i;

		for(i=0; i<nStripe; i++)
		{
			if(locElem[row][i]<0) {
				locElem[row][i] = col;
				elem[row][i] = val;
				return 0;
			}
		}
        return -1;
	}

    //Multiply a scalar value
    void Multiply(real mul)
    {
        int i, j;

        for(i=0; i<dim; i++) {
            for(j=0; j<nStripe; j++) elem[i][j] *= mul;
        }
    }

    //Subtract a matrix
    void Subtract(CSparceMatrix &Mat);

    //Tranpose Mat and Store it
    void Transpose(CSparceMatrix &Mat);

	//Print the matrix
	void Print(ostream &outs)
	{
		int i,j,k;
		int loc;

		outs.precision(5);
		
		//outs.setf(ios::scientific, ios::floatfield);
		outs.setf(ios::fixed, ios::floatfield);
		for(i=0; i<dim; i++) {
			for(j=0; j<dim; j++) {
				loc = -1;
				for(k=0; k<nStripe; k++) 
					if(locElem[i][k]==j) loc = k;
				if (loc >= 0) outs << setw(14) << elem[i][loc] <<"("<<j+1<<")"<< ",";


				//if (loc >= 0) outs << setw(14) << elem[i][loc];// << ",";
				//else       outs << setw(14) << 0.;// << ",";
			}
			outs<<endl;
		}
	}
};

#endif