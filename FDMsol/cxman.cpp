#include <iostream>
#include <string>
#include <stdio.h>

#include "cxman.h"

#include "define.h"
//Constructor of CXTable
CCXTable::CCXTable(void)
{
    int i;

    //Initialize member variables
    for(i=0; i<NUM_GRP; i++) {
        SigSca[i] = 0.;
        SigAbs[i] = 0.;
        SigFis[i] = 0.;
        nuSigFis[i] = 0.;

        SigTot[i] = 0.;
    }
}

//Read CX Table
void CCXTable::ReadCXTable(istream &ins)
{
    int i, j;
    char oneChr, card[CARDLEN];
    bool bEnd = false;

    //Read '('
    ins >> oneChr;

    while(bEnd==false)
    {
        ins >> card;
		if (!strcmp(card, "DiffCoeff")) { //Read nu * fission xs
			for (i = 0; i<NUM_GRP; i++) ins >> DiffCoeff[i];
		}
        else if(!strcmp(card, "SigAbs")) { //Read capture xs
            for(i=0; i<NUM_GRP; i++) ins >> SigAbs[i];
        }
        //else if(!strcmp(card, "SigFis")) { //Read fis. xs
        //    for(i=0; i<NUM_GRP; i++) ins >> SigFis[i];
        //}
        else if(!strcmp(card, "nuSigFis")) { //Read nu * fission xs
            for(i=0; i<NUM_GRP; i++) ins >> nuSigFis[i];
        }
		
		else if (!strcmp(card, "SigChi")) { //Read nu * fission xs
			for (i = 0; i<NUM_GRP; i++) ins >> SigChi[i];
		}
		else if (!strcmp(card, "SigSca")) { //Read scattering xs
			for (i = 0; i<NUM_GRP; i++) {
				SigSca[i] = 0.;
				for (j = 0; j<NUM_GRP; j++) {
					ins >> SigScaDiff[i][j];
					SigSca[i] += SigScaDiff[i][j];
				}
			}
		}
        else if(!strcmp(card, ENDSTR)) bEnd = true;
    }
}

//Read CX library
void CCXManager::ReadCXTables(istream &ins)
{
    int i, id;
    char oneChr, card[CARDLEN];
    bool bEnd = false;

    //Read '('
    ins >> oneChr;

    while(bEnd==false)
    {
        ins >> card;

        if(!strcmp(card, "CXTableNum"))
        {
            ins >> nCXTBL;
            CXTBL = new CCXTable [nCXTBL];
        }
        else if(!strcmp(card, "CXTable")) //Read "CXTable" sub-section
        {
            ins >> id; //Read ID of CXTable
            if( (id>=0) && (id<nCXTBL) ) CXTBL[id].ReadCXTable(ins); // Read "CXTable"
        }
        else if(!strcmp(card, ENDSTR)) bEnd = true;
    }

    //Manipulate CXs
    for(i=0; i<nCXTBL; i++) CXTBL[i].ManipulateCX();
}

//Manipulate CX data
void CCXTable::ManipulateCX(void)
{
    int i;

    for(i=0; i<NUM_GRP; i++) {
        SigTot[i] = SigAbs[i] + SigSca[i];
		SigRmv[i] = SigTot[i] - SigScaDiff[i][i];
    }

}

