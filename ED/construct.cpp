#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<vector>
#include "define.h"
#include <bitset>
#include <sstream>
#include <cstdint>
#include <omp.h> 
int binaryToDecimal(std::string &String);



void constStates(std::vector<std::bitset<Nbit>> &validStates, uint64_t NST) {
    uint64_t gSector = 0;
    
    validStates.clear();
    std::cout << "Starting while loop with NST = " << NST << std::endl;

    #pragma omp parallel
    {
        std::vector<std::bitset<Nbit>> localValidStates;
        uint64_t localGSector = 0;

        #pragma omp for schedule(dynamic) nowait
        for (uint64_t i = 0; i < NST; i++) { 
            std::bitset<Nlink> state(static_cast<unsigned long>(i)); 
            std::bitset<Nbit> newState;

            bool isGauss = true;
            for (int p = 0; p < VOL; p++) {             
                int px = 2 * p, py = 2 * p + 1;
                int nx = 2 * next[DIM - 1][p], ny = 2 * next[DIM - 2][p] + 1;
                int e1 = state[px] ? 1 : -1;
                int e2 = state[py] ? 1 : -1;
                int e3 = state[nx] ? 1 : -1;
                int e4 = state[ny] ? 1 : -1;
                
                int G = (parity[p] == 0) ? 0 : -1;
                int G1 = (parity[p] == 0) ? GL : -GL;
                
                G = -G + (e1 + e2 - e3 - e4) / 2;
                if (GL != 0) G = G1 + G; 

                if (G != 0 && G != 1) {
                    isGauss = false;
                    break;
                }

                newState[3 * p] = G;
                newState[3 * p + 1] = state[px];
                newState[3 * p + 2] = state[py];
            }

            if (isGauss) {
                localValidStates.push_back(newState);
                localGSector++;
            }
        }

        #pragma omp critical
        {
            validStates.insert(validStates.end(), localValidStates.begin(), localValidStates.end());
            gSector += localGSector;
        }
    }

    std::cout << "# of states satisfying Gauss law with occupation either 0 or 1: " << gSector << std::endl;
}

int binaryToDecimal(std::string &String) {
    int decimalValue = 0;
    int base = 1;
    int i;
    for(i=String.size() - 1;i>=0;i--){
        if(String[i] == '1'){
          decimalValue += base;
        }
        base *= 2;
    }
    return decimalValue;
}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------
int checkGL(std::bitset<Nbit> state, int site){
  
  int e1,e2,e3,e4;
  int px, py, nx, ny, p;
  int G;
  p = site;
  px = 3*p+1; py = 3*p+2;
  nx = 3*next[DIM-1][p]+1; ny = 3*next[DIM-2][p]+2;
  e1 = state[px] ? 1 : -1;
  e2 = state[py] ? 1 : -1;
  e3 = state[nx] ? 1 : -1;
  e4 = state[ny] ? 1 : -1;
            
  
  // Put G = 0. find corrsponding fermion occupation # = G .
  if(parity[site] == 0) G = 0;
  else                  G =-1;
  
  G = state[3*p] + G - 0.5*(e1+e2-e3-e4); 
  return G;
}

void compileStates(std::vector<int> &validStates, std::vector<std::bitset<N>> &stateStrings, std::vector<std::bitset<Nbit>> &States){
   int i, p;
   long long int k;
   
   for(i=0;i<validStates.size();i++){
       k = validStates[i];
       std::bitset<Nbit> newState;
       std::bitset<N> StateLinks = stateStrings[i];
       for(p=0;p<VOL;p++){
        newState[3*p] =  conf[k][p];
        newState[3*p+1] = StateLinks[2*p];
        newState[3*p+2] = StateLinks[2*p+1];
       }
       //std::cout << "NewState[" << i << "] = " << newState << std::endl;
       States.push_back(newState);
   }
   
}


*/

 
/* -----------------------------------------------------------------------------------------------------------------------------------------------------------
void generateStates(bool *states){

  long long int i;
  int j,p;
  int val;
  long long int count=0;
  for(i=0;i<NST;i++){
      for(j=0;j<N;j++){
          //std::cout<< i << " " << j << std::endl;
          states[j+N*i] = (i >> j) & 1;
          //std::cout<<states[j+N*i]<< std::endl;
       }

       std::string s;
       for(j=0;j<N;j++){
           s += states[i*N+j]?'1':'0';
       }
       val = binaryToDecimal(s);
       for(j=0;j<N;j++){
          //site index 
          p=j/2;
          if(j%2==0){ xlink[p][val]= states[j+N*i];}
          if(j%2==1){ ylink[p][val]= states[j+N*i];}
       } 
       count++;
  }
  
  std::cout<<"total # of states : " << count<<std::endl;

}

std::vector<std::string> printStates(bool *states) {
    long long int i;
    int j;
    
    int val;
    //size_t l;
    std::vector<std::string> str(NST);
    for(i=0;i<NST;i++){
        std::string s;
        for(j=0;j<N;j++){
            s += states[i*N+j]?'1':'0';
        }
        val = binaryToDecimal(s);
        if(val<NST){
          str[val] = s;
        } 
    }
    return str;

}

void constStates(std::vector<int> &validStates, std::vector<std::string> &stateStrings){
    int p;
    long long int i = 0;
    long long int gSector=0;
    validStates.clear();
    while(i<NST){
        bool isGauss = true;
        for(p=0;p<VOL;p++){ */
            /* Check : If fermion occupation number is either 0 or 1 */
            /*if(checkGL(i, p) != 0 && checkGL(i,p) != 1){
                isGauss = false;
                break;
            }
        }
        if(isGauss){
           validStates.push_back(i);
           gSector++;
        }
        //int decVal = binaryToDecimal(stateStrings[i]);
        for(p=0;p<VOL;p++){
            conf[i][p] = checkGL(i, p);
        }
        i++;
    }
    std::cout<<"# of states satisfying Gauss law with occupation either 0 or 1: " << gSector<<std::endl;
    
    std::vector<std::string> validStateStrings;
    for(const auto &stateIndex : validStates){
        validStateStrings.push_back(stateStrings[stateIndex]);
    }
    stateStrings = validStateStrings;
}

*/

