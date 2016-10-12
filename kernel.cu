#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include "CrsMatrix.h"
#include <cassert>
#include <stdexcept>
#include <vector>
#include <complex>
#include <cmath> 

typedef float VType;

__global__
void kronk_kernel(int* A_Col,int* B_Col,VType* A_Val,VType* B_Val,int* C_Col,VType* C_Val, int na, int counter, int size_kk, int start_k, int start_kk ) {  

  //** Global Mem is slower than local Mem, so local variables are created so Global Var on accessed once
 int s_size_kk ;
 int s_start_kk ;
 int s_start_k ;
 int s_counter;
 int s_na ;

 s_size_kk = size_kk;
 s_start_kk = start_kk;
 s_start_k = start_k;
 s_counter = counter;
 s_na = na;


    int inx = threadIdx.x + blockDim.x * blockIdx.x;
    int iny = threadIdx.y + blockDim.y * blockIdx.y;



     C_Col[inx + iny*s_size_kk + s_counter] = A_Col[iny + s_start_k] + B_Col[inx + s_start_kk] * s_na;
     C_Val[inx + iny*s_size_kk + s_counter] = A_Val[iny + s_start_k] * B_Val[inx + s_start_kk] ;

 }


void CRS_EXTRN_PROD( CrsMatrix<VType>  &C,CrsMatrix<VType> const &A,CrsMatrix<VType> const &B,VType **pd,int **pi) {

   assert(A.row()==A.col());
   assert(B.row()==B.col());
   int n=A.row()*B.row();
   C.resize(n,n);
   int na = A.row();
   int coa = A.nonZero();
   int cob = B.nonZero();
   int i,alpha,beta,counter=0;
   C.resizecv(coa*cob);

   int *d_A_Col,*d_B_Col,*d_C_Col;
   VType *d_A_Val,*d_B_Val,*d_C_Val;

 // Allocating the arrays on the Device
    cudaMalloc((void **) &d_A_Col, sizeof(int) * coa ); 
    cudaMalloc((void **) &d_A_Val, sizeof(VType) * coa );  
    cudaMalloc((void **) &d_B_Col, sizeof(int) * cob ); 
    cudaMalloc((void **) &d_B_Val, sizeof(VType) * cob );  
    cudaMalloc((void **) &d_C_Col, sizeof(int) * coa*cob ); 
    cudaMalloc((void **) &d_C_Val, sizeof(VType) * coa*cob );  

 // Copying the necessary arrays from the Host to the Device
    cudaMemcpy( d_A_Col, &A.colind_[0], sizeof(int) * coa, cudaMemcpyHostToDevice );
    cudaMemcpy( d_A_Val, &A.values_[0], sizeof(VType) * coa, cudaMemcpyHostToDevice );
    cudaMemcpy( d_B_Col, &B.colind_[0], sizeof(int) * cob, cudaMemcpyHostToDevice );
    cudaMemcpy( d_B_Val, &B.values_[0], sizeof(VType) * cob, cudaMemcpyHostToDevice );

  
  for (i=0;i<n;i++) {
   C.setRow(i,counter);
   beta = int(i/na);
   alpha = i - beta * na;

   int size_k = A.getRowPtr(alpha+1) - A.getRowPtr(alpha); 
   int size_kk = B.getRowPtr(beta+1) - B.getRowPtr(beta); 
   int start_k = A.getRowPtr(alpha);
   int start_kk = B.getRowPtr(beta);

   //*** Creating the grid and block geometry and launching the Kernel on the Device
   dim3 grid_dim(  1 ,  1);
   dim3 block_dim(size_kk  , size_k);
   kronk_kernel<<<grid_dim, block_dim>>>(d_A_Col, d_B_Col, d_A_Val, d_B_Val, d_C_Col, d_C_Val, na, counter,
                                              size_kk, start_k, start_kk);
   counter = counter + size_kk*size_k;
  }

  //** Copying the Arrays for CrsMatrix "C" from the Device to the Host
   cudaMemcpy( &C.colind_[0], d_C_Col, sizeof(int) * coa*cob, cudaMemcpyDeviceToHost );
   cudaMemcpy( &C.values_[0], d_C_Val, sizeof(VType) * coa*cob, cudaMemcpyDeviceToHost );
   C.setRow(n,counter);

   pi[0] = d_B_Col;  // 
   pd[0] = d_B_Val;  // Saves these pointers to devie memory so that
   pi[1] = d_C_Col; // they can be used on the next iterations of CRS_EXTRN_PROD_PARTIAL
   pd[1] = d_C_Val; // 

  cudaFree(d_A_Col); // Freeing the Device memory for
  cudaFree(d_A_Val); // the CrsMatrix "A" because its no longer needed

}



void CRS_EXTRN_PROD_PARTIAL( CrsMatrix<VType>  &C,CrsMatrix<VType> const &A,CrsMatrix<VType> const &B,VType **pd,int **pi) {

   assert(A.row()==A.col());
   assert(B.row()==B.col());
   int n=A.row()*B.row();
   C.resize(n,n);
   int na = A.row();
   int coa = A.nonZero();
   int cob = B.nonZero();
   int i,alpha,beta,counter=0;
   C.resizecv(coa*cob);

   int *d_A_Col,*d_B_Col,*d_C_Col;
   VType *d_A_Val,*d_B_Val,*d_C_Val;

   d_B_Col= pi[0];
   d_B_Val= pd[0];
   d_A_Col= pi[1];
   d_A_Val= pd[1];

    cudaMalloc((void **) &d_C_Col, sizeof(int) * coa*cob ); 
    cudaMalloc((void **) &d_C_Val, sizeof(VType) * coa*cob );  

  
  for (i=0;i<n;i++) {
    C.setRow(i,counter);
    beta = int(i/na);
    alpha = i - beta * na;

    int size_k = A.getRowPtr(alpha+1) - A.getRowPtr(alpha); 
    int size_kk = B.getRowPtr(beta+1) - B.getRowPtr(beta); 
    int start_k = A.getRowPtr(alpha);
    int start_kk = B.getRowPtr(beta);

   //*** Creating the grid and block geometry and launching the Kernel on the Device
    dim3 grid_dim(  1 ,  1);
    dim3 block_dim( size_kk , size_k);
    kronk_kernel<<<grid_dim, block_dim>>>(d_A_Col, d_B_Col, d_A_Val, d_B_Val, d_C_Col, d_C_Val, na, counter,
                                            size_kk, start_k, start_kk);
    counter = counter + size_kk*size_k;
  }

  //** Copying the Arrays for CrsMatrix "C" from the Device to the Host
   cudaMemcpy( &C.colind_[0], d_C_Col, sizeof(int) * coa*cob, cudaMemcpyDeviceToHost );
   cudaMemcpy( &C.values_[0], d_C_Val, sizeof(VType) * coa*cob, cudaMemcpyDeviceToHost );
   C.setRow(n,counter);

   pi[1] = d_C_Col; // Saves these pointers to devie memory so that
   pd[1] = d_C_Val;  // they can be used on the next iterations of CRS_EXTRN_PROD_PARTIAL

   cudaFree(d_A_Col); // Freeing the Device memory for
   cudaFree(d_A_Val); // the CrsMatrix "A" because its no longer neede


}


void Free_GPU_MEM(VType **pd,int **pi) {

 //** Freeing the Device memory  **//
     cudaFree(pi[0]); 
     cudaFree(pi[1]);
     cudaFree(pd[0]);
     cudaFree(pd[1]);
}























