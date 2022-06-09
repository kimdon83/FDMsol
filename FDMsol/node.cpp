#include <iomanip>
#include <string>
#include <math.h>

#include "node.h"

#include "define.h"
//Constructor of Cell
CNode::CNode(void)
{
	int i;

	//Initialize member variables
	for (i = 0; i < NEIGH_NUM; i++) NeighborNode[i] = NULL;
	for (i = 0; i < 3; i++) {
		nMesh[i] = 1; meshSize[i] = 1.;
		meshArea[i] = 1.;
	}
	meshVol = 1.;
	cxtblId = -1;
	CX = NULL;
}

void CNode::SetNetLossMatrix(CSparceMatrix * pMatL, real * MeshSize, int loc[3], int kIn, int jIn, int ** nMeshX, int * nMeshXY, int nMeshXYZ, real * Albedo)
{
	int i, j, k, g,eng;
	int pos = 0;
	int nmesh[3];
	real DN[6][NUM_GRP];
	real DM[6][NUM_GRP];
	real beta;
	real mbeta[3][NUM_GRP];
	real area[3];
	real vol;
	CNode *pNode;
	real alb;


	vol = MeshSize[XDIR] * MeshSize[YDIR] * MeshSize[ZDIR];
	area[XDIR] = MeshSize[YDIR] * MeshSize[ZDIR];
	area[YDIR] = MeshSize[XDIR] * MeshSize[ZDIR];
	area[ZDIR] = MeshSize[XDIR] * MeshSize[YDIR];

	//get number of mesh
	for (i = 0; i < 3; i++) nmesh[i] = GetNumberOfMesh(i);

	for (g = 0; g < NUM_GRP; g++) {

		mbeta[XDIR][g] = CX->GetDiffCoeff(g) / MeshSize[XDIR];
		mbeta[YDIR][g] = CX->GetDiffCoeff(g) / MeshSize[YDIR];
		mbeta[ZDIR][g] = CX->GetDiffCoeff(g) / MeshSize[ZDIR];
	}

	for (i = 0; i < 6; i++) {
		pNode = NeighborNode[i];
		alb = Albedo[i];
		for (g = 0; g < NUM_GRP; g++) {
			if (pNode) {
				beta = pNode->CX->GetDiffCoeff(g) / MeshSize[(int)(i / 2)];
				DN[i][g] = 2 * beta * mbeta[(int)(i / 2)][g] / (beta + mbeta[(int)(i / 2)][g])*area[(int)(i / 2)];
			}
			else {
				if (alb >= 1.e10) DN[i][g] = 2 * mbeta[(int)(i / 2)][g] * area[(int)(i / 2)];
				else {
					beta = Albedo[i] / 2;
					DN[i][g] = 2 * mbeta[(int)(i / 2)][g] * beta / (mbeta[(int)(i / 2)][g] + beta)*area[(int)(i / 2)];
				}
			}
		}
	}


	for (k = 0; k < nmesh[2]; k++) {
		for (j = 0; j < nmesh[1]; j++) {
			for (i = 0; i < nmesh[0]; i++) {
				pos = k*nMeshXY[kIn] + j*nMeshX[kIn][jIn] + i + loc[0] + loc[1] + loc[2];

				for (g = 0; g < NUM_GRP; g++) {
					if (i == 0)				DM[XL][g] = DN[XL][g]; else DM[XL][g] = mbeta[XDIR][g] * area[XDIR];
					if (i == nmesh[0] - 1)	DM[XR][g] = DN[XR][g]; else DM[XR][g] = mbeta[XDIR][g] * area[XDIR];
					if (j == 0)				DM[YL][g] = DN[YL][g]; else DM[YL][g] = mbeta[YDIR][g] * area[YDIR];
					if (j == nmesh[1] - 1)	DM[YR][g] = DN[YR][g]; else DM[YR][g] = mbeta[YDIR][g] * area[YDIR];
					if (k == 0)				DM[ZL][g] = DN[ZL][g]; else DM[ZL][g] = mbeta[ZDIR][g] * area[ZDIR];
					if (k == nmesh[2] - 1)	DM[ZR][g] = DN[ZR][g]; else DM[ZR][g] = mbeta[ZDIR][g] * area[ZDIR];

					if ((k == 0) && (NeighborNode[ZL])) {
						pMatL->SetElement(pos + nMeshXYZ * g, pos + nMeshXYZ*g - nMeshXY[kIn - 1],
							-DM[ZL][g]);
					}
					if (k > 0) {
						pMatL->SetElement(pos + nMeshXYZ * g, pos + nMeshXYZ*g - nMeshXY[kIn],
							-DM[ZL][g]);
					}
					if ((j == 0) && (NeighborNode[YL])) {
						pMatL->SetElement(pos + nMeshXYZ * g, pos + nMeshXYZ*g - nMeshX[kIn][jIn - 1],
							-DM[YL][g]);
					}
					if (j > 0) {
						pMatL->SetElement(pos + nMeshXYZ * g, pos + nMeshXYZ*g - nMeshX[kIn][jIn],
							-DM[YL][g]);
					}
					if ((i > 0) || ((i == 0) && (NeighborNode[XL]))) {
						pMatL->SetElement(pos + nMeshXYZ * g, (pos - 1) + nMeshXYZ*g,
							-DM[XL][g]);
					}

					pMatL->SetElement(pos + nMeshXYZ * g, pos + nMeshXYZ * g,
						DM[XL][g] + DM[XR][g] + DM[YL][g] + DM[YR][g] + DM[ZL][g] + DM[ZR][g] + CX->GetCX_RMV(g)*vol);


					if ((i < nMesh[0] - 1) || ((i == nMesh[0] - 1) && (NeighborNode[XR]))) {
						pMatL->SetElement(pos + nMeshXYZ * g, (pos + 1) + nMeshXYZ*g,
							-DM[XR][g]);
					}

					if ((j < nMesh[1] - 1) || ((j == nMesh[1] - 1) && (NeighborNode[YR]))) {
						pMatL->SetElement(pos + nMeshXYZ * g, (pos + nMeshX[kIn][jIn]) + nMeshXYZ*g,
							-DM[YR][g]);
					}


					if ((k < nMesh[2] - 1) || ((k == nMesh[2] - 1) && (NeighborNode[ZR]))) {
						pMatL->SetElement(pos + nMeshXYZ * g, pos + nMeshXY[kIn] + nMeshXYZ * g,
							-DM[ZR][g]);
					}


					//setting scattering matrix
					/*for (eng = 0; eng < NUM_GRP; eng++) {
						if (g != eng) {
							pMatL->SetElement(pos + nMeshXYZ * g, pos + nMeshXYZ * (eng),
								-CX->GetScaDiff(eng, g) *vol);
						}
					}*/
					for (eng = g + 1; eng < NUM_GRP; eng++) {
						pMatL->SetElement(pos + nMeshXYZ * (g + 1), pos + nMeshXYZ * g,
							-CX->GetScaDiff(g, eng) *vol);
					}
				}

			}
		}
	}
}

void CNode::SetFissProdMatrix(CSparceMatrix * pMatF, real * MeshSize, int loc[3], int kIn, int jIn, int ** nMeshX, int * nMeshXY, int nMeshXYZ)
{
	int i, j, k, g;
	int pos = 0;
	real vol;
	int nmesh[3];

	vol = MeshSize[XDIR] * MeshSize[YDIR] * MeshSize[ZDIR];
	//get number of mesh
	for (i = 0; i < 3; i++) nmesh[i] = GetNumberOfMesh(i);

	for (k = 0; k < nmesh[2]; k++) {
		for (j = 0; j < nmesh[1]; j++) {
			for (i = 0; i < nmesh[0]; i++) {
				//set mid position at Matrix
				pos = k*nMeshXY[kIn] + j*nMeshX[kIn][jIn] + i + loc[0] + loc[1] + loc[2];
				for (g = 0; g < NUM_GRP; g++) {
						//Set a diagonal term
						pMatF->SetElement(pos,pos + nMeshXYZ * g,CX->GetCX_NUFIS(g)*vol);
					}
				}
			}
		}
	}
void CNode::SetXiMatrix(CSparceMatrix & Xi, real * MeshSize, int loc[3], int nMeshXTOT, int nMeshXYTOT, int nMeshXYZTOT)
{
	int i, j, k, g;
	int pos = 0;
	real vol;
	int nmesh[3];

	vol = MeshSize[XDIR] * MeshSize[YDIR] * MeshSize[ZDIR];
	//get number of mesh
	for (i = 0; i < 3; i++) nmesh[i] = GetNumberOfMesh(i);

	for (k = 0; k < nmesh[2]; k++) {
		for (j = 0; j < nmesh[1]; j++) {
			for (i = 0; i < nmesh[0]; i++) {
				//set mid position at Matrix
				pos = k*nMeshXYTOT + j*nMeshXTOT + i + loc[0] + loc[1] + loc[2];
				for (g = 0; g < NUM_GRP; g++) {
					Xi.SetElement(pos + nMeshXYZTOT * g, pos ,CX->GetCX_CHI(g));
				}
			}
		}
	}
}


//Constructor of the CStructure class
CStructure::CStructure(void)
{
	//Initialize
	Node = NULL;
	nyNode = 0;
	nxNode = 0;
	xNodeSize = NULL;
	yNodeSize = NULL;
}

//Destructor of the CStructure class
CStructure::~CStructure(void)
{
	int i, j,k;

	//Erase memories
	if (Node) {
		for (k = 0; k < nzNode; k++) {
			for (i = 0; i < nyNode; i++) {
				for (j = 0; j < nxNode; j++) delete Node[k][i][j];
				delete[] Node[k][i];
			}
			delete[] Node[k];
		}
		delete[] Node;
	}
	if (nodeTypeId) {
		for (k = 0; k < nNodeType; k++) {
			delete[] nodeTypeId[k];
		}
		delete[] nodeTypeId;
	}
	if (xNodeSize) delete[] xNodeSize;
	if (yNodeSize) delete[] yNodeSize;
	if (zNodeSize) delete[] zNodeSize;

}

//Read "Geometry" section in a input file
void CStructure::ReadGeometry(istream &ins, int nDim)
{
	int i, j,k;
	char oneChr, card[CARDLEN];
	bool bEnd = false;

	//Read '('
	ins >> oneChr;

	while (bEnd == false)
	{
		ins >> card;

		if (!strcmp(card, "NodeNum"))
		{
			if (nDim == 1) {
				nyNode = 1;
				yNodeSize = new real[nyNode];
				yNodeSize[0] = 1.;
			}
			else if (nDim == 2) {

				ins >> nxNode;
				ins >> nyNode;
			}
			else {

				ins >> nxNode;
				ins >> nyNode;
				ins >> nzNode;
			}
			Node = new CNode ***[nzNode];
			for (j = 0; j < nzNode; j++) {
				Node[j] = new CNode**[nyNode];
				for (i = 0; i < nyNode; i++) {
					Node[j][i] = new CNode* [nxNode];
				}
			}
			
		}
		else if (!strcmp(card, "xNodeSize")) //Read "Size" card
		{
			xNodeSize = new real[nxNode];
			for (i = 0; i < nxNode; i++) ins >> xNodeSize[i];
		}
		else if (!strcmp(card, "yNodeSize")) //Read "Size" card
		{
			yNodeSize = new real[nyNode];
			for (i = 0; i < nyNode; i++) ins >> yNodeSize[i];
		}
		else if (!strcmp(card, "zNodeSize")) //Read "Size" card
		{
			zNodeSize = new real[nzNode];
			for (i = 0; i < nzNode; i++) ins >> zNodeSize[i];
		}



		else if (!strcmp(card, "NodeType"))
		{
			ins >> nNodeType;
			nodeTypeId = new int*[nNodeType];
			for (i = 0; i < nNodeType; i++) {
				nodeTypeId[i] = new int[nzNode];
			}

			for (i = 0; i < nNodeType; i++)
				for (j = 0; j < nzNode; j++)
				{
					ins >> card;
					nodeTypeId[i][j] = atoi(card);
				}
		}
		else if (!strcmp(card, "Configuration"))
		{
			//
				for (i = 0; i < nyNode; i++)
				{
					for (j = 0; j < nxNode; j++)
					{
						ins >> card;
						for (k = 0; k < nzNode; k++) {
							if (strcmp(card, STR_EMPTY)) {

								Node[k][i][j] = new CNode;
								Node[k][i][j]->cxtblId = nodeTypeId[atoi(card)][k];
								//Node[i][j]->configId = atoi(card);
								//Node[i][j]->cxtblId = atoi(card);
															}
							else {
								Node[k][i][j] = NULL;
							}
						}
					}
				}
			//}
		}
		else if (!strcmp(card, ENDSTR)) bEnd = true; //Read ");"
	}

	//Link neighbor nodes
	LinkNeighborNodes();
}

void CStructure::LinkNeighborNodes(void)
{
	int i, j,k;

	for (k = 0; k < nzNode; k++) {
		for (i = 0; i < nyNode; i++) {
			for (j = 0; j < nxNode; j++) {
				if (Node[k][i][j]) {
					//X-LEFT
					if (j > 0)				Node[k][i][j]->LinkNodeXL(Node[k][i][j - 1]);
					//X-RIGHT
					if (j < nxNode - 1)		Node[k][i][j]->LinkNodeXR(Node[k][i][j + 1]);
					//Y-LEFT
					if (i > 0)				Node[k][i][j]->LinkNodeYL(Node[k][i - 1][j]);
					//Y-RIGHT
					if (i < nyNode - 1)		Node[k][i][j]->LinkNodeYR(Node[k][i + 1][j]);
					//Z-LEFT
					if (k > 0)				Node[k][i][j]->LinkNodeZL(Node[k-1][i ][j]);
					//Z-RIGHT
					if (k < nzNode - 1)		Node[k][i][j]->LinkNodeZR(Node[k+1][i ][j]);
				}
			}
		}
	}
}

void CStructure::LinkCXTable(CCXManager *pCXMan)
{
	int i, j, k;

	for (i = 0; i < nyNode; i++) {
		for (j = 0; j < nxNode; j++) {

			//Node[i][j]->CXs = new CCXTable*[nzNode];
			for (k = 0; k < nzNode; k++) {
				if (Node[k][i][j]) {

					Node[k][i][j]->SetCXTable(pCXMan->GetCXTable(Node[k][i][j]->GetCXTableID()));

					//Node[i][j]->CXs[k] = (pCXMan->GetCXTable(nodeTypeId[Node[i][j]->configId][k]));
				}
				//Node[i][j]->SetCXTable(pCXMan->GetCXTable(Node[i][j]->GetCXTableID()));
			//}
			}
		}
	}
}


int CStructure::CalculateNumberOfMeshes(real * MeshSize, int ** nMeshX, int * nMeshXY)
{
	int nx, ny, nz;
	int mnx, mnxy, mnxyz = 0;
	int i, j, k;
	int tot = 0;

	for (k = 0; k < nzNode; k++) {
		mnxy = 0;
		nz = (int)(zNodeSize[k] / MeshSize[ZDIR]);
		for (j = 0; j < nyNode; j++) {
			mnx = 0;
			ny = (int)(yNodeSize[j] / MeshSize[YDIR]);
			for (i = 0; i < nxNode; i++) {
				if (Node[k][j][i]) {
					nx = (int)(xNodeSize[i] / MeshSize[XDIR]);
					mnx += nx;
					Node[k][j][i]->SetMeshNumbers(nx, ny, nz);
				}
			}
			nMeshX[k][j] = mnx;
			mnxy += mnx*ny;
		}
		nMeshXY[k] = mnxy;
		mnxyz += mnxy*nz;
	}

	return mnxyz;
}

void CStructure::SetNetLossMatrix(CSparceMatrix * pMatL, real * MeshSize, int ** nMeshX, int * nMeshXY, int nMeshAll, real * albedo)
{
	int i, j,k ;
	int ny = 0, nz = 0;
	int loc[3] = { 0,0,0 };

	for (k = 0; k < nzNode; k++) {
		for (j = 0; j < nyNode; j++) {
			for (i = 0; i < nxNode; i++) {
				if (Node[k][j][i]) {
					Node[k][j][i]->SetNetLossMatrix(pMatL, MeshSize, loc, k, j, nMeshX, nMeshXY, nMeshAll, albedo);
					loc[0] += Node[k][j][i]->GetNumberOfMesh(XDIR);
				}
				if (i == 0) {
					ny = Node[k][j][i]->GetNumberOfMesh(YDIR);
					nz = Node[k][j][i]->GetNumberOfMesh(ZDIR);
				}
			}
			loc[0] = 0;
			loc[1] += nMeshX[k][j] * ny;
		}
		loc[1] = 0;
		loc[2] += nMeshXY[k] * nz;
	}
}

void CStructure::SetFissProdMatrix(CSparceMatrix * pMatF, real * MeshSize, int ** nMeshX, int * nMeshXY, int nMeshAll)
{
	int i, j,k ;
	int ny = 0, nz = 0;
	int loc[3] = { 0,0,0 };

	for (k = 0; k < nzNode; k++) {
		for (j = 0; j < nyNode; j++) {
			for (i = 0; i < nxNode; i++) {
				if (Node[k][j][i]) {
					Node[k][j][i]->SetFissProdMatrix(pMatF, MeshSize, loc, k, j, nMeshX, nMeshXY, nMeshAll);
					loc[0] += Node[k][j][i]->GetNumberOfMesh(XDIR);
				}
				if (i == 0) {
					ny = Node[k][j][i]->GetNumberOfMesh(YDIR);
					nz = Node[k][j][i]->GetNumberOfMesh(ZDIR);
				}
			}
			loc[0] = 0;
			loc[1] += nMeshX[k][j] * ny;
		}
		loc[1] = 0;
		loc[2] += nMeshXY[k] * nz;
	}
}

void CStructure::SetXiMatrix(CSparceMatrix & Xi, real * MeshSize, int ** nMeshX, int * nMeshXY, int nMeshTot)
{
	int i, j, k;
	int ny = 0, nz = 0;
	int loc[3] = { 0,0,0 };

	for (k = 0; k < nzNode; k++) {
		for (j = 0; j < nyNode; j++) {
			for (i = 0; i < nxNode; i++) {
				if (Node[k][j][i]) {
					Node[k][j][i]->SetXiMatrix(Xi, MeshSize, loc, nMeshX[k][j], nMeshXY[k], nMeshTot);
					loc[0] += Node[k][j][i]->GetNumberOfMesh(XDIR);
				}
				if (i == 0) {
					ny = Node[k][j][i]->GetNumberOfMesh(YDIR);
					nz = Node[k][j][i]->GetNumberOfMesh(ZDIR);
				}
			}
			loc[0] = 0;
			loc[1] += nMeshX[k][j] * ny;
		}
		loc[1] = 0;
		loc[2] += nMeshXY[k] * nz;
	}
}
