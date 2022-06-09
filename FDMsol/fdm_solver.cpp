#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include "fdm_solver.h"
#include "sparce_matrix.h"

#include "define.h"
//Constructor of a FDMSolver object
CFDMSolver::CFDMSolver(void)
{
    //Initialize member variables
    nDim = 1;
    meshSize[XDIR] = 1.;
    meshSize[YDIR] = 1.;
    meshSize[ZDIR] = 1.;
	memset(albedo, 0, sizeof(real)*6);
    ke = 1.e10;
    omega = 1.0;
	method = 0;
    nHigherOrder = 0;

    keff = 1.;
    Flux = NULL;
}

//Destructor of a FDMSolver object
CFDMSolver::~CFDMSolver(void)
{
    //Erase memories
}

//Read an input file
void CFDMSolver::ReadInput(char *filename)
{
    char card[CARDLEN];

    ifstream inFp(filename); //Make input file stream

    //Read input
    while(inFp.peek()!=EOF)
    {
        inFp>>card;

        //Read "Options" section
          if(!strcmp(card, "Options")) this->ReadOptions(inFp);
        //Read "Geometry" section
        else if(!strcmp(card, "Geometry")) Str.ReadGeometry(inFp, nDim);
        //Read "CXLibrary" section
        else if(!strcmp(card, "CXLibrary")) CXMan.ReadCXTables(inFp);
    }

	//Prepare calculation
	Str.LinkCXTable(&CXMan);
}

//Read the calculation options
void CFDMSolver::ReadOptions(istream &ins)
{
    int i;
    char oneChr, card[CARDLEN];
    bool bEnd = false;

    //Read '('
    ins >> oneChr;

    while(bEnd==false)
    {
        ins >> card;

        if(!strcmp(card, "Dimension")) {
            ins >> nDim;
        }
        else if(!strcmp(card,"MeshSize")) {
            for(i=0; i<nDim; i++) ins >> meshSize[i]; //XDIR, YDIR, ZDIR
        }
		else if(!strcmp(card, "Albedo")) {
			for(i=0; i<nDim*2; i++) ins >> albedo[i]; //XL, XR, YL, YR, ZL, ZR
		}
        else if(!strcmp(card, "ke")) {
            ins >> ke;
        }
        else if(!strcmp(card, "omega")) {
            ins >> omega;
        }
        else if(!strcmp(card, "HigherOrder")) {
            ins >> nHigherOrder;
        }
		else if (!strcmp(card, "Method")) {
			ins >> card;
			//if (!strcmp(card, "Jacobi")) method = 0;
			if (!strcmp(card, "GS")) method=1;
			else if (!strcmp(card, "SOR")) method=2;

		}
        else if(!strcmp(card, ENDSTR)) bEnd = true;
    }
}

//Analyze the core problem by the FDM method
void CFDMSolver::AnalyzeCore(void)
{
	int i;
	int nIt;
	int nMeshAll;
	int nStripe;

	int **nMeshX = NULL, *nMeshXY = NULL;
	//allocate memory for number of mesh
	nMeshXY = new int[Str.nzNode];
	memset(nMeshXY, 0, sizeof(int)*Str.nyNode);


	nMeshX = new int*[Str.nzNode];
	for (i = 0; i < Str.nzNode; i++) {
		nMeshX[i] = new int[Str.nyNode]; 
		memset(nMeshX[i], 0, sizeof(int)*Str.nyNode);
	}
	
	//-------------------------------------------------------------------------
	// Allocate memories
	//-------------------------------------------------------------------------
	//Calculation the dimension of the FDM matrix
	nMeshAll=Str.CalculateNumberOfMeshes(meshSize, nMeshX, nMeshXY);
		
	//Set the number of stripes of the net-loss matrix
	nStripe = (1 + nDim*2)+NUM_GRP;
	//Allocate memories of the vector and the matrices
	Flux.SetDimension(nMeshAll*NUM_GRP);
	MatL.SetDimension(nMeshAll*NUM_GRP, nStripe);
	MatF.SetDimension(nMeshAll*NUM_GRP, 2*NUM_GRP-1);
	Xi.SetDimension(nMeshAll*NUM_GRP, 2 * NUM_GRP - 1);
	//-------------------------------------------------------------------------
	// Set matrices
	//-------------------------------------------------------------------------
    SetMatrices(MatL, MatF, meshSize, nMeshX, nMeshXY, nMeshAll);

	//-------------------------------------------------------------------------
	// Solve the eigenvalue equation:
    //     Calculate the fundamental-mode eigenvalue and eigenvector
    //-------------------------------------------------------------------------
    nIt = ExecutePowerIterationMethod(MatL, MatF, Flux, keff);
    cout<<"Number of Iterations = "<<setw(4)<<nIt<<", Keff = "<<keff<<endl;
   
    //Erase memories
	if (nMeshX) {
		for (i = 0; i < Str.nzNode; i++)
			delete[] nMeshX[i];
		delete[] nMeshX;
	}
	if (nMeshXY)
		delete[] nMeshXY;
}

void CFDMSolver::SetMatrices(CSparceMatrix & MatL, CSparceMatrix & MatF, real * MeshSize, int ** nMeshX, int * nMeshXY, int nMeshAll)
{
	Str.SetNetLossMatrix(&MatL, meshSize, nMeshX, nMeshXY, nMeshAll, albedo);
	Str.SetFissProdMatrix(&MatF, meshSize, nMeshX, nMeshXY, nMeshAll);
	Str.SetXiMatrix(Xi, MeshSize, nMeshX, nMeshXY, nMeshAll);
	//Apply the Wielandt method
	if (ke<1.e10) {
		//1/ke * F
		MatF.Multiply(1. / ke);
		//L-1/ke*F
		MatL.Subtract(MatF);

		//Return back to the original F matrix
		MatF.Multiply(ke);

	}
}

int CFDMSolver::ExecutePowerIterationMethod(CSparceMatrix &MatL, CSparceMatrix &MatF,
                                             CVector &Flux, real &keff)
{
	CVector prefisSource;
	CVector fisSource;
	CVector Source; //Source vector
    CVector PrevFlux; //Flux at the previous iteration

    bool bConv = false;
    real prevKeff;
    real maxErr;
    int nIt=-1;
	real norm, denom;

	keff = 1.0;
	//-------------------------------------------------------------------------
	// Allocate memories
	//-------------------------------------------------------------------------
    prefisSource.SetDimension(Flux.GetDimension());
	fisSource.SetDimension(Flux.GetDimension());
	Source.SetDimension(Flux.GetDimension());
	PrevFlux.SetDimension(Flux.GetDimension());

	//-------------------------------------------------------------------------
	// Solve the eigenvalue equation
	//-------------------------------------------------------------------------

	fisSource.Multiply(MatF, Flux);   
	Source.Multiply(Xi,fisSource);

    nIt = 0;
    do
    {
		prefisSource.SetValues(Source);
		
		Source.Multiply(1. / keff);
		prevKeff = keff;

        PrevFlux.SetValues(Flux);

		if (method == 0)	Flux.SolveByPointJacobiIterativeMethod(MatL, Source);
        else if (method == 1) Flux.SolveByGaussSeidelIterativeMethod(MatL, Source);
		else if (method == 2) Flux.SolveBySuccessiveOverRelaxationMethod(MatL, Source, omega);

        fisSource.Multiply(MatF, Flux);
		Source.Multiply(Xi, fisSource);

		denom = fisSource.Multiply(prefisSource);
		norm  = fisSource.Multiply(fisSource);

		keff = prevKeff*norm / denom;

        //if(fabs(prevKeff-keff)<1.e-5) bConv = true;
        maxErr = Flux.GetMaxRelDiff(PrevFlux);
        if(maxErr<1.e-6) bConv = true;
        
        nIt++;
		
		cout << " "<<"Outer Iter : " << nIt << " " << " keff :" << "\t" << keff  ;
		cout << " Error = " << maxErr<<endl;
    }while(!bConv);

	ofstream Flux1("Flux.txt");
	Flux.Print(Flux1);
	Flux1.close();

	ofstream Source1("Source.txt");
	Source.Print(Source1);
	Source1.close();

    //Calcualte keff considering the applied ke
    if(ke<1.e10) keff = 1./(1./keff+1./ke);
	
    return nIt;
}
