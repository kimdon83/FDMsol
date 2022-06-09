#include "fdm_solver.h"

#include "define.h"
int main(void)
{
    CFDMSolver FDSol;

	//time start
	clock_t start_point, end_point;
	start_point = clock();
	double caltime;

    //Read an input
//	FDSol.ReadInput("input_1D_1G_hom_test.txt");	  //		1.35887         
//	FDSol.ReadInput("input_1D_1G_het_test.txt");    //	1.21957				
	//FDSol.ReadInput("input_1D_2G_hom_test.txt");      //	0.51886

	//FDSol.ReadInput("input_IAEA1.txt"); 			//  1.02913
	//FDSol.ReadInput("input_IAEA2.txt");			//  1.02864
	//FDSol.ReadInput("input_IAEA3.txt");			//  1.02887

	FDSol.ReadInput("input.txt");			//  
//	FDSol.ReadInput("input_candu_noDepConDev.txt");			

    //Solve the core analysis problem
//	omp_set_num_threads(4);
	FDSol.AnalyzeCore();

	end_point = clock();
	caltime = (end_point - start_point)/(double)( CLOCKS_PER_SEC);

	cout << caltime << "(sec)" << endl;

	system("pause");
    return 0;
}
