 #include<stdio.h>
 #include<math.h>
 #include<string.h>
 #include<stdlib.h>
 #include "ranlxd.h"
 #include "mkl.h"
 #include "define.h"


void initconf(int FLAG){
  char rfname[100];
  extern void readConf();

  extern int readrand(const char *);

  extern int convert(int, int, int);
  extern void addCluster(CLUSTER *, CLUSTER *, int *, int, int *, int, int);
  snprintf(rfname, sizeof(rfname), "rng_state.bin");

  int i,j,k,pty, ID;
  int ix,iy,it;
  int *sitelist;
  int *linklist;
  CLUSTER *icluster;
  CLUSTER *fcluster;
  sitelist = (int *)mkl_malloc(LT * sizeof(int), alignment);
  totalclusters = 0;
  for(i=0; i<VOL2; i++) breakup[i] = 0;
   /* the initial configuration shown below corresponds to maximally
    * modifiable configuration for different G -sector G=(even,odd)
    * up denotes the orientation of the gauge spins, and (1,0) the
    * staggered fermion occupation number
    */
   i=0;
   if(FLAG == 0){  // G = (0,0)
   for(it=0;it<LT;it++){
     for(iy=0;iy<LY;iy++){
     for(ix=0;ix<LX;ix++){
       pty=(ix+iy)%2;
       if(pty){ // odd site
	      iocc[i]=1; lspin[i]=1;
        if(LY>1) lspin[i+LX*LY*LT]=-1; //vertical link
       }
       else{       // even site
        iocc[i]=0; lspin[i]= -1;
        if(LY>1) lspin[i+LX*LY*LT]=1; //vertical link
       }
       i++;
     }}}   
   }
   else if(FLAG == 1){ // G = (1,-1)
   for(it=0;it<LT;it++){
     for(iy=0;iy<LY;iy++){
     for(ix=0;ix<LX;ix++){
       pty=(ix+iy)%2;
       if(pty){ // odd site
	      iocc[i]=0; lspin[i]=1;
        if(LY>1) lspin[i+LX*LY*LT]=-1; //vertical link
       }
       else{       // even site
        iocc[i]=1; lspin[i]= 1;
        if(LY>1) lspin[i+LX*LY*LT]=-1; //vertical link
       }
       i++;
     }}}
   
   }
   else if(FLAG == 2){  // G = (-1,1)
   for(it=0;it<LT;it++){
     for(iy=0;iy<LY;iy++){
     for(ix=0;ix<LX;ix++){
       pty=(ix+iy)%2;
       if(pty){ // odd site
	      iocc[i]=0; lspin[i]=-1;
        if(LY>1) lspin[i+LX*LY*LT]=-1; //vertical link
       }
       else{       // even site
        iocc[i]=1; lspin[i]= 1;
        if(LY>1) lspin[i+LX*LY*LT]=1; //vertical link
       }
       i++;
     }}}
   
   }
   else if(FLAG == 3){ // G = (2,-2)
   for(it=0;it<LT;it++){
     for(iy=0;iy<LY;iy++){
     for(ix=0;ix<LX;ix++){
       pty=(ix+iy)%2;
       if(pty){ // odd site
	      iocc[i]=1; lspin[i]=1;
        if(LY>1) lspin[i+LX*LY*LT]=1; //vertical link
       }
       else{       // even site
        iocc[i]=0; lspin[i]= -1;
        if(LY>1) lspin[i+LX*LY*LT]=-1; //vertical link
       }
       i++;
     }}}  
   }
   
   else if(FLAG == 4){ // G = (-2,2)
   for(it=0;it<LT;it++){
     for(iy=0;iy<LY;iy++){
     for(ix=0;ix<LX;ix++){
       pty=(ix+iy)%2;
       if(pty){ // odd site
	      iocc[i]=1; lspin[i]=-1;
        if(LY>1) lspin[i+LX*LY*LT]=-1; //vertical link
       }
       else{       // even site
        iocc[i]=0; lspin[i]= 1;
        if(LY>1) lspin[i+LX*LY*LT]=1; //vertical link
       }
       i++;
     }}}
   }
   /** Add a new FLAG for read conf **/
   else if(FLAG == 5){
     readConf();
     readrand(rfname);
   }
   
   
  /*Build the clusters*/
  linklist = NULL;
  icluster = NULL;
  fcluster = NULL;
  for(i=0; i<LX; i++)
  {
    for(j=0; j<LY; j++)
    {
      for(k=0; k<LT; k++)
      {
        sitelist[k] = convert(i,j,k);
      }
      
      addCluster( icluster, fcluster, sitelist, LT, linklist, 0, 0);
      fcluster = firstcluster;
      icluster = firstcluster->prev;
      //printf("cluster site = ");
      /*for(int m = 0; m < LT; m++) {
        
         printf("%d ", sitelist[m]);
      }
      printf("\n"); */
      
     

    }
  }
  mkl_free(sitelist);
  sitelist = NULL;
  linklist = (int *)mkl_malloc(LT * sizeof(int), alignment);
  /*links in the x-direction*/
  for(i=0; i<LX; i++)
  {
    for(j=0; j<LY; j++)
    {
      for(k=0; k<LT; k++)
      {
        linklist[k] = convert(i,j,k);
      }
      addCluster( icluster, fcluster, sitelist, 0, linklist, LT, 0);
      fcluster = firstcluster;
      icluster = firstcluster->prev;
      //printf("cluster link site = ");
      
      /*for(int m = 0; m < LT; m++) {
        
         printf("x direc = %d ", linklist[m]);
      }
      printf("\n");*/
      
    }
  }
  /*links in the y-direction*/
  if(LY > 1)
  {
    for(i=0; i<LX; i++)
    {
      for(j=0; j<LY; j++)
      {
        for(k=0; k<LT; k++)
        {
          linklist[k] = convert(i,j,k)+LX*LY*LT;
        }
        addCluster( icluster, fcluster, sitelist, 0, linklist, LT, 0);
        fcluster = firstcluster;
        icluster = firstcluster->prev;
        //printf("cluster link site = ");
        
        /*for(int m = 0; m < LT; m++) {
        
         printf("y direc = %d ", linklist[m]);
        }
        printf("\n"); */
      }
  }
  }
 
  mkl_free(linklist);
}

void printConf(void)
{
  int i, j, t;

  printf("sites:\n");
  for(t=0; t<LT; t++)
  {
    for(j=0; j<LY; j++)
    {
      for(i=0; i<LX; i++)
      {
        printf("%d ",iocc[i + LX * j + LX * LY * t]);
      }
    }
    printf("\n");
  }
  printf("links horizontal:\n");
  for(t=0; t<LT; t++)
  {
    for(j=0; j<LY; j++)
    {
      for(i=0; i<LX; i++)
      {
        printf("%d ",lspin[i + LX * j + LX * LY * t]);
      }
    }
    printf("\n");
  }
  printf("links vertical:\n");
  for(t=0; t<LT; t++)
  {
    for(j=0; j<LY; j++)
    {
      for(i=0; i<LX; i++)
      {
        printf("%d ",lspin[i + LX * j + LX * LY * t + LX*LY*LT]);
      }
    }
    printf("\n");
  }
}

