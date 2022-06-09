#ifndef _CELL_H_
#define _CELL_H_

#include "define.h"
#include "cxman.h"
#include "sparce_matrix.h"

/**************************************************************************************************
 * Unit Node
**************************************************************************************************/
class CNode {
private:
    CNode *NeighborNode[NEIGH_NUM]; //[NEIGHBOR_INDEX]
	int nMesh[3]; //XDIR, YDIR, ZDIR
	real meshSize[3];
	real meshVol;
    real meshArea[3];

public:
    int cxtblId; //xs table id given in the input file
	int configId; // 
    CCXTable *CX; // This Cell's cross-sections
	//CCXTable **CXs; // adresses of the cross sections of Node

public:
    CNode(void); //Constructor
	~CNode(void) { delete[] CX; delete[] NeighborNode; } //Destructor

	void SetMeshSize(real size[3]) {
		memcpy(meshSize, size, sizeof(real)*3); 
		meshVol = meshSize[XDIR]*meshSize[YDIR]*meshSize[ZDIR];
        meshArea[XDIR] = meshSize[YDIR]*meshSize[ZDIR];
        meshArea[YDIR] = meshSize[XDIR]*meshSize[ZDIR];
        meshArea[ZDIR] = meshSize[XDIR]*meshSize[YDIR];
	}
	real GetMeshSize(int dir) { return meshSize[dir]; }

	real GetMeshVol(void) { return meshVol; }

	void SetMeshNumbers(int nX, int nY, int nZ) { 
		nMesh[XDIR] = nX; nMesh[YDIR] = nY; nMesh[ZDIR] = nZ;
	}
	int GetNumberOfMesh(int dir) { return nMesh[dir]; }
	int GetNumberOfMeshes(void) { return nMesh[XDIR]*nMesh[YDIR]*nMesh[ZDIR]; }

	int GetCXTableID(void) { return cxtblId; }
	void SetCXTable(CCXTable *pCX) { CX = pCX; }
	void SetCXTables(CCXTable *pCX) {  }

	void LinkNodeXL(CNode *pNode) { NeighborNode[XL] = pNode; }
	void LinkNodeXR(CNode *pNode) { NeighborNode[XR] = pNode; }
	void LinkNodeYL(CNode *pNode) { NeighborNode[YL] = pNode; }
	void LinkNodeYR(CNode *pNode) { NeighborNode[YR] = pNode; }
	void LinkNodeZL(CNode *pNode) { NeighborNode[ZL] = pNode; }
	void LinkNodeZR(CNode *pNode) { NeighborNode[ZR] = pNode; }

public:
	void SetNetLossMatrix(CSparceMatrix *pMatL, real *MeshSize, int loc[3],int k, int j, int **nMeshX, int *nMeshXY, int nMeshXYZ, real *Albedo);
	void SetFissProdMatrix(CSparceMatrix *pMatF, real *MeshSize, int loc[3],int k, int j, int **nMeshX, int *nMeshXY, int nMeshXYZ);
	void SetXiMatrix(CSparceMatrix &Xi, real * MeshSize, int loc[3], int nMeshXTOT, int nMeshXYTOT, int nMeshXYZTOT);
};

/**************************************************************************************************
 * Structure 
**************************************************************************************************/
class CStructure {
public:
    CNode ****Node; //[Z][Y][X]
    int nyNode;
    int nxNode;
	int nzNode;

    real *xNodeSize;
    real *yNodeSize;
	real *zNodeSize;

	int **nodeTypeId;
	int nNodeType;

public:
    CStructure(void); //Constructor
    ~CStructure(void); //Destructor
	

	//Link CXTable into each Node
	void LinkNeighborNodes(void);
	void LinkCXTable(CCXManager *pCXMan);

    CNode* GetNode(int zId,int yId, int xId) { return Node[zId][yId][xId]; }

public:
    void ReadGeometry(istream &ins, int nDim); //Read "Geometry" section in a input file

	//Functions for making matrices
	//Return the dimension of the system matrix
	int CalculateNumberOfMeshes(real *MeshSize, int **nMeshX, int *nMeshXY);
	//Set a net-loss matrix
	void SetNetLossMatrix(CSparceMatrix *pMatL, real *MeshSize, int **nMeshX
								, int *nMeshXY, int nMeshAll, real *albedo);
	//Set a fission-production matrix
	void SetFissProdMatrix(CSparceMatrix *pMatF, real *MeshSize, int **nMeshX
		, int *nMeshXY, int nMeshAll);
	void SetXiMatrix(CSparceMatrix &ML, real *MeshSize, int **nMeshX, int *nMeshXY, int nMeshTot);

};
#endif