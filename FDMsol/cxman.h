#ifndef _CX_MANAGER_H_
#define _CX_MANAGER_H_

#include "define.h"

/**************************************************************************************************
 * CX Table 
**************************************************************************************************/
class CCXTable {
//-------------------------------------------------------------------------------------------------
// CX Data from Library File
//-------------------------------------------------------------------------------------------------
private:
    real SigSca[NUM_GRP], SigScaDiff[NUM_GRP][NUM_GRP]; // scattering xs
    real SigAbs[NUM_GRP]; // absorption xs
    real SigFis[NUM_GRP]; // fission xs
    real nuSigFis[NUM_GRP]; // nu * ( fission xs )
	real DiffCoeff[NUM_GRP]; //diffusion coefficient
	real SigChi[NUM_GRP];

//-------------------------------------------------------------------------------------------------
// Manipulated CX Data
//-------------------------------------------------------------------------------------------------
private:
    real SigTot[NUM_GRP]; // Total CX
	real SigRmv[NUM_GRP]; // Removal CX

public:
    CCXTable(void); //Constructor
    ~CCXTable(void) {} //Destructor

	real GetCX_CHI(int engGrp) { return SigChi[engGrp]; }
    real GetCX_TOT(int engGrp) { return SigTot[engGrp]; }
    real GetCX_NUFIS(int engGrp) { return nuSigFis[engGrp]; }
    real GetCX_FIS(int engGrp) { return SigFis[engGrp]; }
    real GetCX_ABS(int engGrp) { return SigAbs[engGrp]; }
	real GetCX_RMV(int engGrp) { return SigRmv[engGrp]; }
	real GetScaDiff(int engGrp, int engGrp2) { return SigScaDiff[engGrp][engGrp2]; }
    real GetNu(int engGrp) { return nuSigFis[engGrp] / SigFis[engGrp]; }
	real GetDiffCoeff(int engGrp) { return DiffCoeff[engGrp]; }

public:
    void ReadCXTable(istream &ins);
    void ManipulateCX(void);
};

/**************************************************************************************************
 * CX Manager
**************************************************************************************************/
class CCXManager {
private:
    int nCXTBL;
    CCXTable *CXTBL;

public:
    CCXManager(void) { nCXTBL = 0; CXTBL = NULL; } //Constructor
    ~CCXManager(void) { delete [] CXTBL; } //Destructor

    CCXTable* GetCXTable(int id) { return &CXTBL[id]; }

public:
    void ReadCXTables(istream &ins); //Read "CXTables" section in a input file
};

#endif