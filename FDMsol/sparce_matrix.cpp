#include "sparce_matrix.h"

#include "define.h"

void CSparceMatrix::Subtract(CSparceMatrix &Mat)
{
    int i, j, k;
    bool bMatch;

    //Check dimensions
	if( dim!=Mat.GetDimension() ) {
		cout<<"Dimension mismatch at CSparceMatrix::Subtract"<<endl;
		exit(0);
	}

    //Perform *this - Mat
    for(i=0; i<dim; i++)
    {
        for(j=0; j<Mat.GetNumberOfStripes(); j++)
        {
            bMatch = false;
            for(k=0; k<this->nStripe; k++)
            {
                if(Mat.locElem[i][j]==this->locElem[i][k]) {
                    this->elem[i][k] -= Mat.elem[i][j];
                    bMatch =true;
                    break;
                }
            }
            if(!bMatch) {
                cout<<"Location mismatch at CSparceMatrix::Subtract"<<endl;
                exit(0);
            }
        }
    }
}

void CSparceMatrix::Transpose(CSparceMatrix &Mat)
{
    int i, j;

    //Check dimensions
    if( (dim!=Mat.GetDimension()) || (nStripe!=Mat.GetNumberOfStripes()) ) {
		cout<<"Dimension or stripe number mismatch at CSparceMatrix::Tranpose"<<endl;
		exit(0);
	}

    //Perform Transpose
    for(i=0; i<Mat.dim; i++)
    {
        for(j=0; j<Mat.nStripe; j++)
        {
            if(Mat.locElem[i][j]>=0) {
                if(SetElement(Mat.locElem[i][j], i, Mat.elem[i][j])<0) {
                    cout<<"Error at CSparceMatrix::Tranpose"<<endl;
		            exit(0);
                }
            }
        }
    }
}