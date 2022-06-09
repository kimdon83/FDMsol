#ifndef _FDM_SOLVER_SHJ_H_
#define _FDM_SOLVER_SHJ_H_

#include "define.h"
#include "cxman.h"
#include "node.h"

#include "vector.h"

/**************************************************************************************************
 * FDM Solver
**************************************************************************************************/
class CFDMSolver {
//-------------------------------------------------------------------------------------------------
// Objects for the core calculations
//-------------------------------------------------------------------------------------------------
public:
    CCXManager CXMan; //CX Manager object
    CStructure Str; //Structure object

	CSparceMatrix MatL; //Net-Loss Matrix
	CSparceMatrix MatF; //Fission Matrix
	CSparceMatrix Xi;
//-------------------------------------------------------------------------------------------------
// Calculation Options
//-------------------------------------------------------------------------------------------------
private:
    int nDim;
    real meshSize[3]; //in the order of XDIR, YDIR, ZDIR
	real albedo[6];
    real ke; //expected keff for the Wieland method
    real omega; //extrapolation parameter for the SOR method

    int nHigherOrder; //the maximum order of eigenvectors to calculate
	int method;

//-------------------------------------------------------------------------------------------------
// Output
//-------------------------------------------------------------------------------------------------
private:
    real keff;
    CVector Flux;

public:
    CFDMSolver(void); //Constructor
    ~CFDMSolver(void); //Destructor

//-------------------------------------------------------------------------------------------------
// Methods
//-------------------------------------------------------------------------------------------------
public:
    void ReadInput(char *filename); //Read Input file
        void ReadOptions(istream &ins); //Read "Options" section in a input file

	void AnalyzeCore(void); //Analyze the core problem by the FDM method
        void SetMatrices(CSparceMatrix &MatL, CSparceMatrix &MatF,
			real *MeshSize, int **nMeshX, int *nMeshXY, int nMeshAll);
        int ExecutePowerIterationMethod(CSparceMatrix &MatL, CSparceMatrix &MatF, CVector &Flux, real &keff);
      /*  void CalculateHigherOrderSolution(CSparceMatrix &MatL,  CSparceMatrix &MatF,  CVector *FwdFlux,
                                          CSparceMatrix &MatLT, CSparceMatrix &MatFT, CVector *AdjFlux,
                                          real *eigenvalue);*/
};

#endif