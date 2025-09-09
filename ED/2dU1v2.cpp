/** Exact diagonalization code for matter coupled with spin-1/2 U(1) gauge links using Lanczos algorithm **/

/** constuction Hamiltonian matrix H_1 = -t(\psi^+_x S^+_xy \psi_y + \psi^+_y S^-_xy \psi_x) + t(n_x - 0.5)(n_y - 0.5) + t s^3_xy (n_y - n_x)
                                        + V (n_x - 0.5)(n_y - 0.5) **/ 

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>
#include "mkl.h"
#include <string>
#include <bitset>
#include <sstream>
#include <fstream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <cstdint> 

int *next[2*DIM+1];
int *parity;
int LX,LY,VOL,GL,iter,FLAG;
double Vn;

//bool *states;
uint64_t NST;
typedef Eigen::SparseMatrix<double> SparseMatrix;

//double *Hmatrix;
int alignment;
double *psiBpsi, *Sx, *Sy, *cdw, *plaqOp1;

int main(){
  FILE *fptr;
  char string[50];
  char Hfname[100], rfname[100];
  int i;

  alignment = 64;
  extern void initneighbor(void);
  extern void neighchk();
  extern void constStates(std::vector<std::bitset<Nbit>> &, uint64_t NST);
  extern void constructH1Sparse(std::vector<std::bitset<Nbit>> &, Eigen::SparseMatrix<double> &);
  extern void printMatrix(std::vector<std::bitset<Nbit>> &, Eigen::SparseMatrix<double> &);
  extern void lanczos(Eigen::SparseMatrix<double> &, int , Eigen::VectorXd &, Eigen::MatrixXd &);
  extern void obs(std::vector<std::bitset<Nbit>> &, double *, Eigen::MatrixXd &, Eigen::VectorXd &, int , int, double *, double *, double *, double *,FILE *);
  extern void measure(std::vector<std::bitset<Nbit>> &, double *, double *, double *, double *, double *);
  extern void blockLanczos(Eigen::SparseMatrix<double> &, int , int , Eigen::VectorXd &, Eigen::MatrixXd &);
  extern void spectralanczos(Eigen::SparseMatrix<double> &, int , Eigen::VectorXd &, Eigen::MatrixXd &);
  extern void writeH(char *, Eigen::SparseMatrix<double> &);
  extern void readH(char *, Eigen::SparseMatrix<double> &);
  
  fptr = fopen("QUEUE2D","r");
  if(fptr == NULL){
      std::cout << "Could not open QUEUE2D file" << std::endl;
      exit(1);
  }
  fscanf(fptr,"%s %d\n",string,&LX);
  fscanf(fptr,"%s %d\n",string,&LY);
  fscanf(fptr,"%s %d\n",string,&GL);
  fscanf(fptr,"%s %le\n",string,&Vn);
  fscanf(fptr,"%s %d\n",string,&iter);
  fscanf(fptr,"%s %d\n",string,&FLAG);
  
  //fscanf(fptr,"%s %le\n",string,&V);
  fclose(fptr);
  
  VOL = LX*LY;
  
  if(GL==1)  printf("Code for the Gauss Law sector: (Ge, Go) = (1,-1)  \n");
  if(GL==2)  printf("Code for the Gauss Law sector: (Ge, Go) = (2,-2)  \n");
  if(GL==-1) printf("Code for the Gauss Law sector: (Ge, Go) = (-1,1)  \n");
  if(GL==0)  printf("Code for the Gauss Law sector: (Ge, Go) = (0,0)  \n");
  
  
  sprintf(Hfname,"HmatrixG%dV%.2lfL%d-%d", GL, Vn, LX, LY);
  
  
  //total Hilbert space from Gauge fields and fermions
  //NST = pow(4,VOL);
  uint64_t NST = pow(4, VOL);

  std::cout << "NST = " << NST << std::endl;

  //printf("NST = %lld\n", NST);
  
  /* Initialize nearest neighbours */
  for(i=0;i<=2*DIM;i++){
    next[i] = (int *)calloc(VOL*sizeof(int), alignment);
  }
  parity = (int *)calloc(VOL*sizeof(int), alignment);
  
    
  
  initneighbor();
  //neighchk();

  
  std::vector<std::bitset<Nbit>> validStates;
  int threadList[] = {2,4,8,12,16,20};
  int nTests = sizeof(threadList)/sizeof(threadList[0]);

  for(int t=0;t<nTests;t++){
       int threads = threadList[t];
       omp_set_num_threads(threads);
       double start = omp_get_wtime();
       constStates(validStates, NST);
       double end = omp_get_wtime();

       printf("Threads: %2d | Time: %f seconds\n", threads, end - start);
  }
  //constStates(validStates, NST);
  
  psiBpsi = (double *)calloc(validStates.size() *sizeof(double), alignment);
  cdw     = (double *)calloc(validStates.size() *sizeof(double), alignment);
  Sx      = (double *)calloc(validStates.size() *sizeof(double), alignment);
  Sy      = (double *)calloc(validStates.size() *sizeof(double), alignment);
  plaqOp1 = (double *)calloc(validStates.size() *sizeof(double), alignment);
  
  measure(validStates, psiBpsi, Sx, Sy, cdw, plaqOp1);
  if(FLAG==0){
     /* Generate Sparse matrix (row, col, val) */
     Eigen::SparseMatrix<double> Hmatrix(validStates.size(), validStates.size());

     constructH1Sparse(validStates, Hmatrix);
     std::cout << "End of construction of Hmatrix (sparse)." << std::endl;
     std::cout << "Non zero elements : " << Hmatrix.nonZeros() << std::endl;
     //printMatrix(validStates, Hmatrix);
       
     writeH(Hfname, Hmatrix);
  }
  
  if(FLAG == 1){
     Eigen::SparseMatrix<double> Hmatrix(validStates.size(), validStates.size());
     sprintf(rfname,"HmatrixG%dV%.2lfL%d-%d", GL, Vn, LX, LY);
     readH(rfname, Hmatrix);
     /* Diagonalize Hmatrix with Lanczos algorithm */
     Eigen::VectorXd evals;
     Eigen::MatrixXd evecs; 
     std::cout << "Lanczos algorithm to diagonalize sparse matrix" << std::endl;
     lanczos(Hmatrix, iter, evals, evecs);
     std::cout<<"diagonalization is done." << std::endl;
 
     char filename[100];
     sprintf(filename, "stateG%dV%.2lf-L%d-%dlattice.dat", GL, Vn, LX, LY);
     fptr = fopen(filename, "w");
     obs(validStates, psiBpsi, evecs, evals, validStates.size(), iter,  Sx, Sy, cdw, plaqOp1, fptr); 
     fclose(fptr);
  } 
  /* Clear memory */
  for(i=0;i<=2*DIM;i++){  free(next[i]); }
  free(parity);
  free(psiBpsi);
  free(Sx);
  free(Sy);
  free(cdw);
  free(plaqOp1);
  
  return 0;
}




