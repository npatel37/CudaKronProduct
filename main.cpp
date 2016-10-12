#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <vector>
#include <stdlib.h>
#include <sys/time.h>
#include <ctime>
#include "CrsMatrix.h"


#define CLOCKS_PER_MS (CLOCKS_PER_SEC / 1000)

typedef float VType;

using namespace std;

void CRS_EXTRN_PROD( CrsMatrix<VType>  &C,CrsMatrix<VType> const &A,CrsMatrix<VType> const &B,VType **pd,int **pi);
void CRS_EXTRN_PROD_PARTIAL( CrsMatrix<VType>  &C,CrsMatrix<VType> const &A,CrsMatrix<VType> const &B,VType **pd,int **pi);
void Free_GPU_MEM(VType **pd,int **pi);

class onesite {
public:
	onesite() {
		initialize();
	}
	CrsMatrix<VType> Splus, Sminus, Sz, HamSys;
        VType *pd[2];
        int *pi[2]; 
	void initialize(){
		Matrix<VType> tmp_Splus, tmp_Sminus, tmp_Sz;
		tmp_Splus.resize(2,2); tmp_Splus(0,1) = 1.0;
		tmp_Sminus.resize(2,2); tmp_Sminus(1,0) = 1.0;
		tmp_Sz.resize(2,2);
		tmp_Sz(0,0) = 0.5;
		tmp_Sz(1,1) = -0.5;

		fullMatrixToCrsMatrix(Splus, tmp_Splus);
		fullMatrixToCrsMatrix(Sminus, tmp_Sminus);
		fullMatrixToCrsMatrix(Sz,tmp_Sz);

		CrsMatrix<VType> twositeHam1, twositeHam2, twositeHam3;
		CRS_EXTRN_PROD(twositeHam1, Sz,Sz,pd,pi);
		CRS_EXTRN_PROD(twositeHam2, Splus,Sminus,pd,pi);
		CRS_EXTRN_PROD(twositeHam3, Sminus,Splus,pd,pi);
		//		twositeHam1.DensePrint();
		//		twositeHam2.DensePrint();
		//		twositeHam3.DensePrint();

		//    HamSys.add(HamSys,twositeHam1,1.0);
		HamSys+=twositeHam1;
		HamSys+=twositeHam2;
		HamSys+=twositeHam3;
		//		HamSys.DensePrint();
	}
};

int main()
{
	//~ clock_t begin = clock();
        //unsigned begin = clock();
	Matrix<VType> tmp_Splus, tmp_Sminus, tmp_Sz;
	CrsMatrix<VType> Splus, Sminus, Sz;
	onesite site;
        VType *pd[2]; // this array holds d_B_Col,d_B_Val,d_A_Col,d_A_Val 
        int *pi[2];   // on entry to CRS_EXTRN_PROD_PARTIAL and returns d_B_Col,d_B_Val,d_C_Col,d_C_Val  
	CrsMatrix<VType> testTmp1,testTmp2;
	CrsMatrix<VType> testTmp3,testTmp4;
	CrsMatrix<VType> testTmp5,testTmp6;
	CrsMatrix<VType> testTmp7,testTmp8;
	CrsMatrix<VType> testTmp9;
	
  //** The Full external product is called first to save the pointers for the next call
	CRS_EXTRN_PROD(testTmp1, site.HamSys,site.HamSys, pd,pi); // 4-sites
 //** The array of pointers p now has the device pointers to avoid unnecessary data transfer 
	CRS_EXTRN_PROD_PARTIAL(testTmp2, testTmp1,site.HamSys, pd,pi); // 6 sites
           testTmp2.DensePrint();
	CRS_EXTRN_PROD_PARTIAL(testTmp3, testTmp2,site.HamSys, pd,pi); // 8 sites
	CRS_EXTRN_PROD_PARTIAL(testTmp4, testTmp3,site.HamSys, pd,pi); // 10 sites
	CRS_EXTRN_PROD_PARTIAL(testTmp5, testTmp4,site.HamSys, pd,pi); // 12 sites
	CRS_EXTRN_PROD_PARTIAL(testTmp6, testTmp5,site.HamSys, pd,pi); // 14 sites
	CRS_EXTRN_PROD_PARTIAL(testTmp7, testTmp6,site.HamSys, pd,pi); // 16 sites
	CRS_EXTRN_PROD_PARTIAL(testTmp8, testTmp7,site.HamSys, pd,pi); // 18 sites
	//unsigned begin = clock();
	CRS_EXTRN_PROD_PARTIAL(testTmp9, testTmp8,site.HamSys, pd,pi); // 20 sites	
        unsigned begin = clock();
        CRS_EXTRN_PROD_PARTIAL(testTmp1, testTmp9,site.HamSys, pd,pi); // 22 sites
 

	//site.HamSys.DensePrint();	
	//~ testTmp = CRS_EXTRN_PROD(testTmp,site.HamSys); // 8 sites
	//~ testTmp = CRS_EXTRN_PROD(testTmp,site.HamSys); // 10 sites
	//~ testTmp = CRS_EXTRN_PROD(testTmp,site.HamSys); // 12 sites
	//~ testTmp = CRS_EXTRN_PROD(testTmp,site.HamSys); // 14 sites
	//~ testTmp = CRS_EXTRN_PROD(testTmp,site.HamSys); // 16 sites
	//~ testTmp = CRS_EXTRN_PROD(testTmp,site.HamSys); // 18 sites
	//~ testTmp = CRS_EXTRN_PROD(testTmp,site.HamSys); // 20 sites

	//testTmp.DensePrint();

	//~ clock_t end = clock();
	//unsigned end = clock() / (CLOCKS_PER_SEC);
	//VType elapsed_secs = VType(end - begin);
	//cout << "elapsed time: " << elapsed_secs << std::endl;
	

	clock_t timeElapsed = clock() - begin;
	unsigned msElapsed = timeElapsed / CLOCKS_PER_MS;
	cout << "elapsed time: " << msElapsed << std::endl;
  
        Free_GPU_MEM(pd,pi);
   
	return 0;
}











