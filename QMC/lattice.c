#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ranlxd.h"
#include "define.h"

void lattice(){
   extern int convert(int,int,int);
   int i,p,ix,iy,it;
   int ixp1,ixm1,iyp1,iym1,itp1,itm1;
   int trot;
   trot = 4;
   //not every link touches a plaquette, so this starts them all out at -1.
   if(LY==1){
      for(i=0; i<VOL; i++) link[i] = -1;
   }
   else{
      for(i=0; i<2*VOL; i++) link[i] = -1;
   }
   i=0; p=0;
   for(it=0;it<LT;it++){
   for(iy=0;iy<LY;iy++){
   for(ix=0;ix<LX;ix++){

     /* serial to lattice-coordinates */
     ixc[i]=ix; itc[i]=it; 
     iyc[i]=iy;
     /* indexing for bond variables */
     //lbc[i]=ix; ltc[i]=it;
     /* staggered factors */
     istag[i]=((ix + iy)%2);
     parity[i] = (ix + iy + it)%2;

     /* neighbours */
     itp1=(it+1)%LT;  itm1=(it-1+LT)%LT;
     ixp1=(ix+1)%LX;  ixm1=(ix-1+LX)%LX;
     if(LY > 1) {
      iyp1=(iy+1)%LY; iym1=(iy-1+LY)%LY;
     }

     /* site neighbours forward and backward in time */
     itup[i]=convert(ix,iy,itp1);  itdn[i]=convert(ix,iy,itm1);

     /* bond neighbors forward and backward in time */
     ltup[i]=convert(ix,iy,itp1);  ltdn[i]=convert(ix,iy,itm1);
     if(LY > 1)//bonds in the y-direction
     {
         ltup[i+LX*LY*LT]=convert(ix,iy,itp1)+LX*LY*LT; 
         ltdn[i+LX*LY*LT]=convert(ix,iy,itm1)+LX*LY*LT;
     }

     /* lattice spatial neighbors */
     xp1[i] = convert(ixp1,iy,it);
     xm1[i] = convert(ixm1,iy,it);
     if(LY > 1)
     {
      yp1[i] = convert(ix,iyp1,it);
      ym1[i] = convert(ix,iym1,it);
     }

     /* spatial bonds in and out  a site, in the direction of forward time */
     if(it%trot==0){
         /* ___   ___
            ___   ___
            ___   ___
            ___   ___
         (0,0)
            lout connects i with the other site in the active placquette
            at the same timeslice
         */
         if(ix%2==0){
		       iout[i]=convert(ixp1,iy,it);   iin[i] =convert(ixm1,iy,it);
		       lout[i]=convert(ix,iy,it);     lin[i] =convert(ixm1,iy,it);
             if(LY > 1)
             {
               iin[i] =convert(ix,iym1,it);
               lin[i] =convert(ix,iym1,it)+LX*LY*LT;
               if(iy%2 ==1)
               {
                  iin[i] =convert(ix,iyp1,it);
                  lin[i] =convert(ix,iy,it)+LX*LY*LT;
               }
             }
	       }
	       if(ix%2==1){
		       iout[i]=convert(ixm1,iy,it);   iin[i] =convert(ixp1,iy,it);
		       lout[i]=convert(ixm1,iy,it);   lin[i] =convert(ix,iy,it);
             if(LY > 1)
             {
               iin[i] =convert(ix,iym1,it);
               lin[i] =convert(ix,iym1,it)+LX*LY*LT;
               if(iy%2 ==1)
               {
                  iin[i] =convert(ix,iyp1,it);
                  lin[i] =convert(ix,iy,it)+LX*LY*LT;
               }
             }
	       }
     }
     if(it%trot==1){
         /*    ___   ___
               ___   ___
               ___   ___
               ___   ___
         (0,0)
            lout connects i with the other site in the active placquette
            at the same timeslice
         */
         if(ix%2==0){
		       iout[i]=convert(ixm1,iy,it);   iin[i] =convert(ixp1,iy,it);
		       lout[i]=convert(ixm1,iy,it);   lin[i] =convert(ix,iy,it);
	       }
	       if(ix%2==1){
		       iout[i]=convert(ixp1,iy,it);   iin[i] =convert(ixm1,iy,it);
		       lout[i]=convert(ix,iy,it);     lin[i] =convert(ixm1,iy,it);
	       }
     }
     //find the sites to the left and the right of a link in the x-direction
     lsitel[convert(ix,iy,it)] = i; lsiter[convert(ix,iy,it)] = convert(ixp1,iy,it);
     /*these will only be called for 2D lattices*/
     if(it%trot==2){
         /*        
            |  |  |  |
                   
            |  |  |  |
         (0,0)
            lout connects i with the other site in the active placquette
            at the same timeslice
         */
         if(iy%2==0){
		       iout[i]=convert(ix,iyp1,it);        
		       lout[i]=convert(ix,iy,it)+LX*LY*LT; //vertical links
             if(ix%2 == 0)
             {
               iin[i] =convert(ixm1,iy,it);
               lin[i] =convert(ixm1,iy,it); 
             }
             else
             {
               iin[i] =convert(ixp1,iy,it);
               lin[i] =convert(ix,iy,it); 
             }
	       }
	       if(iy%2==1){
		       iout[i]=convert(ix,iym1,it); 
		       lout[i]=convert(ix,iym1,it)+LX*LY*LT; //vertical links
             if(ix%2 == 0)
             {
               iin[i] =convert(ixm1,iy,it);
               lin[i] =convert(ixm1,iy,it); 
             }
             else
             {
               iin[i] =convert(ixp1,iy,it);
               lin[i] =convert(ix,iy,it); 
             }
	       }
     }
     if(it%trot==3){
         /* |  |  |  | 
                     
            |  |  |  |
               
         (0,0)
            lout connects i with the other site in the active placquette
            at the same timeslice
         */
         if(iy%2==0){
		       iout[i]=convert(ix,iym1,it);          iin[i] =convert(ix,iyp1,it);
		       lout[i]=convert(ix,iym1,it)+LX*LY*LT; lin[i] =convert(ix,iy,it)+LX*LY*LT; //vertical links
	       }
	       if(iy%2==1){
		       iout[i]=convert(ix,iyp1,it);        iin[i] =convert(ix,iym1,it);
		       lout[i]=convert(ix,iy,it)+LX*LY*LT; lin[i] =convert(ix,iym1,it)+LX*LY*LT; //vertical links
	       }
     }
     if(LY > 1)
     {
      /*the sites to the back and front of a link in the y-direction*/
      lsitel[convert(ix,iy,it)+LX*LY*LT] = i; 
      lsiter[convert(ix,iy,it)+LX*LY*LT] = convert(ix,iyp1,it);
     }
     
     /* label the interaction plaquettes */
     if(it%trot==0){
        if(ix%2==0){
	       	ifwd[i]=convert(ix,iy,it);
            if(LY == 1) ibwd[i]=convert(ixm1,iy,itm1);
            else
            {
               if(iy%2 == 0) ibwd[i]=convert(ix,iym1,itm1);
               else ibwd[i]=convert(ix,iy,itm1);
            }
	      }
        if(ix%2==1){
	       	ifwd[i]=convert(ixm1,iy,it);    
            if(LY==1) ibwd[i]=convert(ix,iy,itm1);
            else
            {
               if(iy%2 == 0) ibwd[i] = convert(ix,iym1,itm1);
               else ibwd[i] = convert(ix,iy,itm1);
            }
	      }
     }
     if(it%trot==1){
        if(ix%2==0){
	       	ibwd[i]=convert(ix,iy,itm1);    ifwd[i]=convert(ixm1,iy,it);
	      }
        if(ix%2==1){
	       	ibwd[i]=convert(ixm1,iy,itm1);  ifwd[i]=convert(ix,iy,it);
	      }
     }
     /*THese will only be called for 2D lattices*/
     if(it%trot==2){
        if(iy%2==0){
	       	ifwd[i]=convert(ix,iy,it);
            if(LY == 1) ibwd[i]=convert(ix,iym1,itm1); //plaquettes parallel to y-axis
            else{
               if(ix%2 == 0) ibwd[i] = convert(ixm1,iy,itm1);
               else ibwd[i] = convert(ix,iy,itm1);
            }
	      }
        if(iy%2==1){
	       	ifwd[i]=convert(ix,iym1,it);    
            if(LY == 1) ibwd[i]=convert(ix,iy,itm1);
            else{
               if(ix%2 == 0) ibwd[i] = convert(ixm1,iy,itm1);
               else ibwd[i] = convert(ix,iy,itm1);
            }
	      }
     }
     if(it%trot==3){
        if(iy%2==0){
	       	ibwd[i]=convert(ix,iy,itm1);    ifwd[i]=convert(ix,iym1,it); //plaquettes parallel to y-axis
	      }
        if(iy%2==1){
	       	ibwd[i]=convert(ix,iym1,itm1);  ifwd[i]=convert(ix,iy,it);
	      }
     }
     /* decide the breakup type on the interaction plaquettes */
     if(((ix+it)%2 == 0) &&((it % trot) != 2) && ((it % trot) != 3)){
	 /*  3 ---- 5 ---- 2
	     |             |
	     |             |
	     0 ---- 4 ---- 1
	 */
         plaq[0][p] = i;
         plaq[1][p] = convert(ixp1,iy,it);
         plaq[2][p] = convert(ixp1,iy,itp1);
         plaq[3][p] = convert(ix,iy,itp1);
         plaq[4][p] = i;
         plaq[5][p] = convert(ix,iy,itp1);

	 /* Plaquette p is the active plaquette directly forward in time (0)
	    and backward in time (1) that touches the site i */
	 /*    plaquette p
	  1,i+hat{t} -----1,i+\hat{x}+\hat{t}
	     |             |
	     |             |
	    0,i --------- 0,i+\hat{x}
	 */
         site[0][i] = p;
         site[0][convert(ixp1,iy,it)] = p;
         site[1][convert(ixp1,iy,itp1)] = p;
         site[1][convert(ix,iy,itp1)] = p;
      /*    plaquette p
	     ----ix, it+1---
	     |             |
	     |             |
	     ---i=(ix,iy)---
	    */
         link[i] = p;
         link[convert(ix,iy,itp1)] = p;
         p++;
     } // bracket closing (ix+it)%2==0
     /*only relevant for 2D lattice*/
      if(((iy+it)%2 == 0) &&(((it % trot) == 2) || ((it % trot) == 3))){
	   /*  3 ---- 5 ---- 2
	     |             |
	     |             |
	     0 ---- 4 ---- 1
	   */
         plaq[0][p] = i;
         plaq[1][p] = convert(ix,iyp1,it);
         plaq[2][p] = convert(ix,iyp1,itp1);
         plaq[3][p] = convert(ix,iy,itp1);
         plaq[4][p] = i+LX*LY*LT; //link in y-direction
         plaq[5][p] = convert(ix,iy,itp1)+LX*LY*LT;
	   /* Plaquette p is the active plaquette directly forward in time (0)
	    and backward in time (1) that touches the site i */
	   /*    plaquette p
	   1,i+hat{t} -----1,i+\hat{x}+\hat{t}
	     |             |
	     |             |
	    0,i --------- 0,i+\hat{x}
	   */
         site[0][i] = p;
         site[0][convert(ix,iyp1,it)] = p;
         site[1][convert(ix,iyp1,itp1)] = p;
         site[1][convert(ix,iy,itp1)] = p;
      /*    plaquette p
	     ----ix, it+1---
	     |             |
	     |             |
	     ---i=(ix,iy)---
	    */
         link[i+LX*LY*LT] = p;
         link[convert(ix,iy,itp1)+LX*LY*LT] = p;
         p++;
     } // bracket closing (ix+it)%2==0
     i++;
   }}}
  
  /*for(p=0; p<VOL; p++)
     {
       printf("site = %d, +up = %d, +dn = %d, iout = %d, iin = %d, fwd = %d, bwd = %d \n", p, itup[p], itdn[p], 
       iout[p], iin[p], ifwd[p], ibwd[p]);
     }
   for(p=0; p<VOL; p++)
     {
       printf("link site = %d, +up = %d, +dn = %d, lout = %d, lin = %d lsitel = %d lsiter = %d \n", p, ltup[p], ltdn[p], 
       lout[p], lin[p], lsitel[p], lsiter[p]);
     
     
       for(i=0;i<6;i++)
         {
          printf("site = %d, active plaquette sites = %d \n", p, plaq[i][p]);     
         } 
     }   */  
   /*for(p=0; p<VOL2; p++)
   {
     printf("plaquette %d is made up of (%d,%d,%d,%d), (%d,%d)\n",p,plaq[0][p],plaq[1][p],plaq[2][p],
                                    plaq[3][p],plaq[4][p],plaq[5][p]);
   }  
    for(p=0;p<VOL;p++){
        printf("p =%d site 0 = %d, site 1 = %d itup = %d itdn =%d iout = %d iin=%d\n", p, site[0][p], site[1][p], itup[p], itdn[p], iout[p], iin[p]);
    
    }*/
    /*for(p=0;p<VOL;p++){
        printf("p =%d link = %d, ltup = %d ltdn =%d lout = %d lin=%d\n", p, link[p], ltup[p], ltdn[p], lout[p], lin[p]);
    
    }*/
 }
