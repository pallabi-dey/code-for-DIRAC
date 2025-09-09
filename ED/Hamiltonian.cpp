#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>
#include "mkl.h"
#include<string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unordered_map>

int binaryToDecimal(std::string &String);

typedef Eigen::SparseMatrix<double> SparseMatrix;

/** constuction Hamiltonian matrix H_1 = -t(\psi^+_x S^+_xy \psi_y + \psi^+_y S^-_xy \psi_x) + t(n_x - 0.5)(n_y - 0.5) + t s^3_xy (n_y - n_x)
                                        + V (n_x - 0.5)(n_y - 0.5) **/ 
void constructH1Sparse(std::vector<std::bitset<Nbit>> &validStates, Eigen::SparseMatrix<double> &Hmatrix){     
   int p;
   size_t i,j, idx;
   long long int newval; 
   
   std::vector<Eigen::Triplet<double>> tripletList;
    
   // Create a hash map for fast state lookup
   std::unordered_map<std::bitset<Nbit>, size_t> stateIndexMap;
   for(idx=0;idx<validStates.size();idx++){
       stateIndexMap[validStates[idx]] = idx;
   }
   /*
   for(idx=0;idx<validStates.size();idx++){
      printf("state %llu with %zu\n", validStates[idx].to_ullong(), stateIndexMap[validStates[idx]]);

   }*/
   for(i = 0; i < validStates.size(); i++){
       std::bitset<Nbit> currentState = validStates[i];
       //std::cout<<"cuurentstate: " << currentState << std :: endl;
       for(p=0;p<VOL;p++){
          if((currentState[3*p] == 0 && currentState[3*next[DIM + 1][p]] == 0)) tripletList.emplace_back(i, i, (Vn * 0.25 - 0.5 * tn));
          if((currentState[3*p] == 0 && currentState[3*next[DIM + 2][p]] == 0)) tripletList.emplace_back(i, i, (Vn * 0.25 - 0.5 * tn));
          if((currentState[3*p] == 1 && currentState[3*next[DIM + 1][p]] == 1)) tripletList.emplace_back(i, i, (Vn * 0.25 - 0.5 * tn));
          if((currentState[3*p] == 1 && currentState[3*next[DIM + 2][p]] == 1)) tripletList.emplace_back(i, i, (Vn * 0.25 - 0.5 * tn));
          
          if((currentState[3*p] == 0 && currentState[3*next[DIM + 1][p]] == 1) && currentState[3*p+1] == 1){
            tripletList.emplace_back(i, i, (-Vn * 0.25 + 0.5 * tn));
          }  
          if((currentState[3*p] == 0 && currentState[3*next[DIM + 2][p]] == 1) && currentState[3*p+2] == 1){
            tripletList.emplace_back(i, i, (-Vn * 0.25 + 0.5 * tn));
          }   
          if((currentState[3*p] == 0 && currentState[3*next[DIM + 1][p]] == 1) && currentState[3*p+1] == 0){
            tripletList.emplace_back(i, i, (-Vn * 0.25 - 0.5 * tn));
          }   
          if((currentState[3*p] == 0 && currentState[3*next[DIM + 2][p]] == 1) && currentState[3*p+2] == 0){
            tripletList.emplace_back(i, i, (-Vn * 0.25 - 0.5 * tn));
          }  
          
          if((currentState[3*p] == 1 && currentState[3*next[DIM + 1][p]] == 0) && currentState[3*p+1] == 1){
             tripletList.emplace_back(i, i, (-Vn * 0.25 - 0.5 * tn));
          }   
          if((currentState[3*p] == 1 && currentState[3*next[DIM + 2][p]] == 0) && currentState[3*p+2] == 1){
             tripletList.emplace_back(i, i, (-Vn * 0.25 - 0.5 * tn));
          }   
          if((currentState[3*p] == 1 && currentState[3*next[DIM + 1][p]] == 0) && currentState[3*p+1] == 0){
             tripletList.emplace_back(i, i, (-Vn * 0.25 + 0.5 * tn));
          }   
          if((currentState[3*p] == 1 && currentState[3*next[DIM + 2][p]] == 0) && currentState[3*p+2] == 0){
             tripletList.emplace_back(i, i, (-Vn * 0.25 + 0.5 * tn));
          }  
          
         
          if(currentState[3*p] == 0 && currentState[3*next[DIM + 1][p]] == 1 && currentState[3*p+1] == 0){
             std::bitset<Nbit> newState = currentState;              
             newState[3*p] = 1;
             newState[3*next[DIM + 1][p]] = 0;
             newState[3*p+1] = 1;
             auto it = stateIndexMap.find(newState);
             if(it != stateIndexMap.end()){
                tripletList.emplace_back(i, it->second, -tn);
             }     
            
          } 
          if(currentState[3*p] == 0 && currentState[3*next[DIM + 2][p]] == 1 && currentState[3*p+2] == 0){              
             std::bitset<Nbit> newState = currentState; 
             newState[3*p] = 1;
             newState[3*next[DIM + 2][p]] = 0;
             newState[3*p+2] = 1;
             auto it = stateIndexMap.find(newState);
             if(it != stateIndexMap.end()){
                tripletList.emplace_back(i, it->second, -tn);
             }
          }
          if(currentState[3*p] == 1 && currentState[3*next[DIM + 1][p]] == 0 && currentState[3*p+1] == 1){ 
             std::bitset<Nbit> newState = currentState; 
             newState[3*p] = 0;
             newState[3*next[DIM + 1][p]] = 1;
             newState[3*p+1] = 0;
             auto it = stateIndexMap.find(newState);
             if(it != stateIndexMap.end()){
                tripletList.emplace_back(i, it->second, -tn);
             }
          }
          if(currentState[3*p] == 1 && currentState[3*next[DIM + 2][p]] == 0 && currentState[3*p+2] == 1){                
             std::bitset<Nbit> newState = currentState;
             newState[3*p] = 0;
             newState[3*next[DIM + 2][p]] = 1;
             newState[3*p+2] = 0;
             auto it = stateIndexMap.find(newState);
             if(it != stateIndexMap.end()){
                tripletList.emplace_back(i, it->second, -tn);
             }
          
          }
          
       }
       //fptr << std::endl; 
          
   }
   Hmatrix.setFromTriplets(tripletList.begin(), tripletList.end());
   Hmatrix.makeCompressed();

}

void printMatrix(std::vector<std::bitset<Nbit>> &validStates, Eigen::SparseMatrix<double> &Hmatrix){
    size_t i,j;
    
    std::cout << "Hamiltonian Matrix:" << std::endl;
    for(i=0;i<Hmatrix.outerSize();i++){
        for(Eigen::SparseMatrix<double>::InnerIterator it(Hmatrix, i); it; ++it){
            std::cout << "H(" << it.row() << ", " << it.col() << ") = " 
                      << std::setprecision(6) << it.value() << std::endl;
        }
    }
}
 

/* Diagonalize the sparse Hmatrix with Lanczos algorithm */

void lanczos(Eigen::SparseMatrix<double> &Hmatrix, int k, Eigen::VectorXd &evals, Eigen::MatrixXd &evecs){
    int i,j,n,m;
    double alpha_j, beta_j;
    n = Hmatrix.rows();
    
    Eigen::VectorXd v0 = Eigen::VectorXd::Zero(n);     // Initialize with zero vector
    Eigen::VectorXd v1 = Eigen::VectorXd::Random(n);   // initial guess for Lanczos randomly
    
    //std::cout << "v0 (initial zero vector):\n" << v0 << "\n\n";
    //std::cout << "v1 (random vector):\n" << v1 << "\n\n";
   
    v1.normalize();
    
    //std::cout << "v1 (normalized random vector):\n" << v1 << "\n\n";
   
    std::vector<double> alpha;                          // Diagonal elements of tridiagonal matrix
    std::vector<double> beta;                           // Off-diagonal elements of tridiagonal matrix

    Eigen::MatrixXd V(n, k);                            // Matrix to store Lanczos basis vectors
    V.col(0) = v1;

    for(j=0;j<k;j++){
        Eigen::VectorXd w = Hmatrix*V.col(j);         
        alpha_j = V.col(j).dot(w);                      // Diagonal element
        alpha.push_back(alpha_j);
        w -= alpha_j * V.col(j);                        // Orthogonalize against previous basis vector

        if(j>0){
           w -= beta[j-1] * V.col(j-1);                
        }
        /* Gram-Schmidt orthonormalization */
        for(m=0;m<=j;m++){
            Eigen::VectorXd vm = V.col(m);
            w -= vm.dot(w) * vm;                       
        }
        
        beta_j = w.norm();                              // Off-diagonal element
        if(j<k-1){
            beta.push_back(beta_j);
            V.col(j+1) = w/beta_j;                      // Normalize and store the next Lanczos vector
        }
    }
    
    /* Construction of the tridiagonal matrix */
   Eigen::MatrixXd T = Eigen::MatrixXd::Zero(k, k);
    for(i=0;i<k;i++){
        T(i, i) = alpha[i];
        if(i<k-1){
            T(i, i+1) = beta[i];
            T(i+1, i) = beta[i];
        }
    }

    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(T);
    if(solver.info() != Eigen::Success){
       std::cerr << "Lanczos eigenvalue decomposition failed!" << std::endl;
       exit(1);
    }

   
    evals = solver.eigenvalues();
    Eigen::MatrixXd T_evecs = solver.eigenvectors();
    // Print eigenvalues
    std::cout << "Eigenvalues at GS:" << evals[0] << std::endl;
    
    evecs = V.leftCols(k) * T_evecs;
}



/*------------------------------------------------------------------------------------------------------------------------------------------------------------

void diagonalizeMatrix(double *Hmatrix, int Hsize, double *eigenvalues, double *evecs){
   
    int i, j;

    std::copy(Hmatrix, Hmatrix+Hsize*Hsize, evecs);

    LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', Hsize, evecs, Hsize, eigenvalues);

    // Print eigenvalues
    std::cout << "Eigenvalues at GS:" << eigenvalues[0] << std::endl;*/
    /*for(i=0;i<Hsize;i++){
        std::cout << eigenvalues[i] << " \n";
    }
    std::cout << std::endl;
    */
    // Print eigenvectors
    /*
    std::ostringstream filename;
    filename << "L2-2/Veq" 
             << (V < 0 ? "m" : "") << std::fixed << std::setprecision(6) << std::abs(V)
             << "t/eigenVector-G" 
             << (GL < 0 ? "m" : "") << std::abs(GL) << ".txt";
    std::ofstream fptr(filename.str());
    fptr << "Eigenvectors (each column is an eigenvector):" << std::endl;
    for(i=0;i<Hsize;i++){
        for(j=0;j<Hsize;j++){
            fptr << evecs[i*Hsize+j] << " ";
        }
        fptr << std::endl;
    }
    */
    

//}





