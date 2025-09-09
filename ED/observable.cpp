#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <iomanip>
#include<vector>
#include "define.h"
#include <fstream>
#include <Eigen/Sparse>
#include <Eigen/Dense>


/** Zero temperature case **/
void obs(std::vector<std::bitset<Nbit>> &validStates, double *psiBpsi, Eigen::MatrixXd &evecs, Eigen::VectorXd &evals, int Hsize, int iter, double *Sx, double *Sy, double *cdw, double *plaqOp1, FILE *fptr){
   int i,j;
   double denum     = 0.0;
   double evSq      = 0.0;
   double chiralCnd = 0.0;
   
   double Ex  = 0.0;
   double Ey  = 0.0;
   double CDW = 0.0;   
   double PlaqOp1 = 0.0;
   
   std::cout << std::fixed << std::setprecision(20); 
   fprintf(fptr, "Eval            psiBpsi          Ex          Ey         CDW      flip_plaq\n");
   for(i=0;i<iter;i++){
       evSq = 0.0;
       chiralCnd = 0.0;
       Ex = Ey = 0.0; CDW = 0.0;
       PlaqOp1 = 0.0;
       double norm = 0.0;
       for(j=0;j<Hsize;j++){
           double prob = evecs(j, i) * evecs(j, i);
           evSq      += prob;
           chiralCnd += psiBpsi[j]*prob;
           Ex        += Sx[j]*prob;
           Ey        += Sy[j]*prob;
           CDW       += cdw[j]*prob;
           PlaqOp1   += plaqOp1[j]*prob;
       }
       norm = evSq; 
       chiralCnd /= norm;
       Ex /= norm;
       Ey /= norm;
       CDW /= norm;
       PlaqOp1 /= norm;
       fprintf(fptr,"%le %le %le %le %le %le\n", evals[i], chiralCnd, Ex, Ey, CDW, PlaqOp1);
       denum += evSq;
       if(i==0) printf("zero temp: <psi^+ psi>_s =  %le\n\n", chiralCnd);

   }
   denum /= Hsize;
   std::cout << "denum = " << denum << std::endl;
   //std::cout << "<psi^+ psi> = " << occN/denum << std::endl;
}



void measure(std::vector<std::bitset<Nbit>> &validStates, double *psiBpsi, double *Sx, double *Sy, double *cdw, double *plaqOp1){
  int i,p;
  for(i = 0; i < validStates.size(); i++){
       std::bitset<Nbit> currentState = validStates[i];
       /* store chiral condensate and staggered gauge links */
       psiBpsi[i] = Sx[i] = Sy[i] = 0.0; plaqOp1[i] = 0.0;
       for(p=0;p<VOL;p++){
           int px = next[DIM+1][p]; int py = next[DIM+2][p];
           psiBpsi[i] += (1-2*parity[p])*currentState[3*p];
           Sx[i]   += (currentState[3*p+1] == 0) ? -0.5 : 0.5;
           Sy[i]   += (currentState[3*p+2] == 0) ? -0.5 : 0.5; 
           cdw[i]  += (1-2*parity[p])*(currentState[3*p]-0.5)*((currentState[3*px]-0.5)+(currentState[3*py]-0.5));
           
           double e1 = (currentState[3*p+1] == 0) ? -0.5 : 0.5;
           double e2 = (currentState[3*px+2] == 0) ? -0.5 : 0.5;
           double e3 = (currentState[3*py+1] == 0) ? -0.5 : 0.5;
           double e4 = (currentState[3*p+2] == 0) ? -0.5 : 0.5;
           plaqOp1[i] += (((0.5+e1)*(0.5+e2)*(0.5-e3)*(0.5-e4))+((0.5-e1)*(0.5-e2)*(0.5+e3)*(0.5+e4)));
          
       }
       psiBpsi[i] /= ((double)VOL);
       Sx[i]   /= ((double)VOL);
       Sy[i]   /= ((double)VOL);
       cdw[i]  /= ((double)VOL);
       plaqOp1[i] /= ((double)VOL);
       
  }
}  
