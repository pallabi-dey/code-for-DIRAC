#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"

void initneighbor(void)
{
  int x, y, p;
  extern int convert(int, int);
  /* lattice linear co-ordinates */
  for(y=0;y<LY;y++){
      for(x=0;x<LX;x++){
         p = convert(x,y);
         parity[p] = ((x+y)%2);
         next[DIM][p] = convert(x,y);
         next[DIM+1][p] = convert(((x+1)%LX), y);
         next[DIM+2][p] = convert(x, ((y+1)%LY));
         next[DIM-1][p] = convert(((x-1+LX)%LX), y);
         next[DIM-2][p] = convert(x, ((y-1+LY)%LY));
      
      }
  }
  
  
}

int convert(int x, int y){  
    return x+y*LX;
}

void neighchk() {
    int p;
    for (p = 0; p < VOL; ++p) {
        std::cout << "site = " << p
                  << ", +x = " << next[DIM + 1][p]
                  << ", +y = " << next[DIM + 2][p]
                  << ", -x = " << next[DIM - 1][p]
                  << ", -y = " << next[DIM - 2][p] 
                  << ", parity = " << parity[p] << std::endl;
    }
}

