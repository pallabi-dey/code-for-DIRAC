/* Code to simulate the (2+1)-d fermions coupled to U(1) gauge fields
   multi cluster version 
*/

 #include<stdio.h>
 #include<math.h>
 #include<string.h>
 #include<stdlib.h>
 #include "ranlxd.h"
 #include "mkl.h"
 #include "define.h"


/* Parameters                      *
 * (LX,LT)    : lattice dimension  *
 * J          : coupling constant  *
 * beta       : inverse temperature*
 * eps        : Trotter size       *
 * ieq        : equilibrium iter   *
 * imeas      : # of measurements  *
 * iskp       : # of iter/meas     */


int main(int argc, char *argv[])
{
  if(argc != 2)
    {
      printf("usage %s jobid\n",argv[0]);
      exit(1);
    }
  extern void rlxd_init(int, int);
  extern void lattice();
  extern void neighchk();
  extern void initconf(int);
  extern void printConf(void);
  extern void printconf();
  extern void updateBreakup();
  extern void Cluster(int *, int);
  extern void printClusters(void);
  extern void updateConfig(void);
  extern void resetID(void);
  extern void checkconf();
  extern void improvedMeasure();
  extern void storeSector(void);
  extern double GLsquared(int);
  extern void measure(void);
  extern void fprintConf(FILE *);
  extern void configObs(FILE *);
  extern void readConf();
  extern void writeConf();
  extern int writerand(const char *);
  extern int readrand(const char *);

  int i, j, ieq, imeas;
  char string[60];
  int readval;
  int ntot;
  char st[100];
  char rfname[100], rfname1[100];
  FILE *fptr, *fptrC, *fptr0, *fptr1, *fptr2;
  
  alignment = 64;
  sscanf(argv[1],"%d",&arg);
  // monitor the expectation value of the fermion
  // occupation number and the Sz gauge fields


  /* read file */
  snprintf(string, sizeof(char)*60, "QUEUE1d-%d", arg);
  //snprintf(rfname, sizeof(rfname), "rng_state.bin");
  
  fptr = fopen(string, "r");
  if (fptr == NULL)
  {
    printf("QUEUE error.\n");
    exit(1);
  }
  readval = fscanf(fptr, "%s %d\n", st, &LX);
  if (readval == -1)
    printf("Error\n");
  readval = fscanf(fptr, "%s %d\n", st, &LY);
  if (readval == -1)
    printf("Error\n");
  readval = fscanf(fptr, "%s %d\n", st, &LT);
  if (readval == -1)
    printf("Error\n");
  readval = fscanf(fptr, "%s %lf\n", st, &beta);
  if (readval == -1)
    printf("Error\n");
  readval = fscanf(fptr, "%s %lf\n", st, &t);
  if (readval == -1)
    printf("Error\n");
  readval = fscanf(fptr, "%s %d\n", st, &ieq);
  if (readval == -1)
    printf("Error\n");
  readval = fscanf(fptr, "%s %d\n", st, &imeas);
  if (readval == -1)
    printf("Error\n");
  /* choose the Hamiltonian */
  readval = fscanf(fptr, "%s %d\n", st, &FLAG);
  if (readval == -1)
    printf("Error\n");
  readval = fscanf(fptr, "%s %d\n", st, &histFlag);
  if (readval == -1)
    printf("Error\n");
  readval = fscanf(fptr, "%s %d\n", st, &ETERMS);
  if (readval == -1)
    printf("Error\n");
  readval = fscanf(fptr, "%s %d\n", st, &CONFIG);
  if (readval == -1)
    printf("Error\n");
  readval = fscanf(fptr, "%s %d\n", st, &INDEX);
  if (readval == -1)
    printf("Error\n");
  readval = fscanf(fptr, "%s %d\n", st, &SEED);
  if (readval == -1)
    printf("Error\n");

  fclose(fptr);

  eps = beta / LT;
  if(LY==1) LT = 2 * LT;
  else LT = 4 * LT;
  VOL = LX * LY * LT;
  VOL2 = VOL / 2;
  V2 = 2 * VOL + 1;
  VOL4 = VOL / 4;
  V4 = 4 * VOL + 1;
  SVOL = LX*LY;

  // prob is A / (A + D) (note expression is the same whether E-terms there or not)
  double tx = 1.0;
  double ty = t;
  prob1 = 2.0 / (exp(2.0 * eps * tx) + 1.0);
  prob2 = 2.0 / (exp(2.0 * eps * ty) + 1.0);
  //prob  = 2.0 / (exp(2.0 * eps * (sqrt(tx*ty)) + 1.0));
  printf("2-d spin-1/2 Gauge-Fermion Model. Multi-cluster version.\n");
  printf("LX=%d, LY=%d, 4*LT=%d, beta=%.3f\n", LX, LY, LT, beta);
  printf("Probability to switch interactions = %lf (X) and %lf (Y) \n", prob1, prob2);
  
  snprintf(rfname, sizeof(rfname), "rng_state.bin");
  snprintf(rfname1, sizeof(rfname1), "rng_state-%d.bin", INDEX);


  /* initialize ranlux */
  if(FLAG == 5) readrand(rfname);
  else          rlxd_init(1, SEED);

  /* allocate memory */
  ixc = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  iyc = (int *)mkl_malloc(VOL * sizeof(int), alignment); 
  itc = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  iocc = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  
  if(LY==1) lspin = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  else      lspin = (int *)mkl_malloc(2*VOL * sizeof(int), alignment);
  
  itup = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  itdn = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  xp1 = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  xm1 = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  yp1 = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  ym1 = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  parity = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  iin = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  iout = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  ifwd = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  ibwd = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  lout = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  lin = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  if(LY==1)
  {
    lclus = (int *)mkl_malloc(VOL * sizeof(int), alignment);
    ltup = (int *)mkl_malloc(VOL * sizeof(int), alignment);
    ltdn = (int *)mkl_malloc(VOL * sizeof(int), alignment);
    lsitel = (int *)mkl_malloc(VOL * sizeof(int), alignment);
    lsiter = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  }
  else
  {
    lclus = (int *)mkl_malloc(2*VOL * sizeof(int), alignment);
    ltup = (int *)mkl_malloc(2*VOL * sizeof(int), alignment);
    ltdn = (int *)mkl_malloc(2*VOL * sizeof(int), alignment);
    lsitel = (int *)mkl_malloc(2*VOL * sizeof(int), alignment);
    lsiter = (int *)mkl_malloc(2*VOL * sizeof(int), alignment);
  }
  iclus = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  istag = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  cflag = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  breakup = (int *)mkl_malloc(VOL2 * sizeof(int), alignment);
  
  if(LY==1)
    link = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  else
    link = (int *)mkl_malloc(2 * VOL * sizeof(int), alignment);
  
  for (i = 0; i < 6; i++)
  {
    plaq[i] = (int *)mkl_malloc(VOL2 * sizeof(int), alignment);
  }
  for (i = 0; i < 2; i++)
    site[i] = (int *)mkl_malloc(VOL * sizeof(int), alignment);
  
  
  int *list;
  int size = 0;
  idmax = 0;
  list = NULL;


  lattice();

  initconf(FLAG);
  checkconf();
  /*printConf();
  printconf();
  printClusters();*/
  //updateBreakup();
  /*extern int plaquettecategory(int);
  for(int p=0;p<VOL2;p++){
     printf("category = %d\n", plaquettecategory(p));
  
  }*/
  //printClusters();
  for(i=0;i<ieq;i++){ 
    //printf("equli:%d\n", i);
    //printConf();
    //if(i%1000) printconf();
    checkconf();
    //printconf();
    updateBreakup();
    //printClusters();
    updateConfig();
    resetID();
  }
  printf("equli done.\n");
  if(FLAG == 5){
    readConf();
    //readrand(rfname);
       
  }
  snprintf(string, sizeof(char)*60, "OnePoint-beta%.1fty%.2f.dat", beta, ty);
  fptr0 = fopen(string,"w");
  //snprintf(string, sizeof(char)*60, "Configs-with-Gsq16-beta%.1f.dat", beta);
  //fptr1 = fopen(string,"w");
  if(SVOL<4){
     snprintf(string, sizeof(char)*60, "eflux-fermN-profile-beta%.1fL-%d-%d.dat", beta, LX, LY);
     fptr2 = fopen(string,"w");
     fprintf(fptr2, " Site\tOcc         xlink       ylink        plaqObs\n");
  }
  for(i=0;i<=imeas;i++){ 
    checkconf();
    //printconf();
    //printConf();
    updateBreakup();
    //printClusters();
    updateConfig();
    resetID();
    measure();
    storeSector();
    G2 = GLsquared(i);
    fprintf(fptr0,"%le %le %le %le %le %le\n", psiBpsi, eFlux, G2, plaqOp1, stagEfluxX, stagEfluxY); 
    /*if(sector2 == 1 ) fprintf(fptr1,"======================Config = %d with Gsq = %0.1lf ==============\n",i, G2);
    if(sector2 == 1 ) fprintConf(fptr1);*/
    
    if(SVOL<4){
      configObs(fptr2); 
      //printConf(); 
    
    }
    if(i%1000 == 0){
      writeConf();
      writerand(rfname1);
    }  
  }
  fclose(fptr0);
  //fclose(fptr1);
  
  /* free variables */
  mkl_free(ixc);   mkl_free(itup); mkl_free(ltup); 
  mkl_free(iyc);   mkl_free(itdn); mkl_free(ltdn); 
  mkl_free(itc);   mkl_free(iout); mkl_free(ifwd);
  mkl_free(lclus); mkl_free(iin);  mkl_free(ibwd);
  mkl_free(iclus); mkl_free(lin);  mkl_free(lout);
  mkl_free(iocc);  mkl_free(lspin); 
  mkl_free(lsiter);
  mkl_free(lsitel);
  
  
  mkl_free(xp1); mkl_free(yp1);
  mkl_free(xm1); mkl_free(ym1);

  mkl_free(parity);
  mkl_free(istag);
  mkl_free(cflag);

  mkl_free(breakup);
  
  
  mkl_free(link);
  for(i=0;i<6;i++) mkl_free(plaq[i]);
  for(i=0;i<2;i++) mkl_free(site[i]);
  return 0;
}

int convert(int x,int y,int t){
   int n;
   n = t*LX*LY + y*LX + x;
   return n;
}



void checkconf(){

  int i, p;
  /*for(p=0; p<VOL2; p++)
   {
     printf("plaquette %d is made up of (%d,%d,%d,%d)=(%d,%d,%d,%d), (%d,%d)=(%d,%d)\n",p,plaq[0][p],plaq[1][p],plaq[2][p],
                                    plaq[3][p], iocc[plaq[0][p]],iocc[plaq[1][p]],
       iocc[plaq[2][p]],iocc[plaq[3][p]],  plaq[4][p],plaq[5][p], lspin[plaq[4][p]],lspin[plaq[5][p]]);
   }
   */
  for(p=0;p<VOL2;p++){
      if(iocc[plaq[0][p]] == iocc[plaq[1][p]] && iocc[plaq[1][p]] == iocc[plaq[2][p]] && iocc[plaq[2][p]] == iocc[plaq[3][p]]){
         if(lspin[plaq[4][p]] != lspin[plaq[5][p]]){
          printf("error: invalid plaquette = %d, type0!\n", p);
          exit(1);
         }      
      }
     else if(iocc[plaq[0][p]] != iocc[plaq[3][p]]){
       if(lspin[plaq[4][p]] != lspin[plaq[5][p]] && iocc[plaq[0][p]] == iocc[plaq[2][p]] && iocc[plaq[0][p]] != iocc[plaq[1][p]]){
          if((iocc[plaq[0][p]] == 1 && lspin[plaq[4][p]] == 1) || ((iocc[plaq[0][p]] == 0 && lspin[plaq[4][p]] == 0))) continue;
       }
       else{ printf("error : invalid plaquette = %d, type1;\n", p); exit(1);}         
     } 
      
       else if(iocc[plaq[0][p]] == iocc[plaq[3][p]] && iocc[plaq[1][p]] == iocc[plaq[2][p]] && iocc[plaq[0][p]] != iocc[plaq[1][p]]){
        if((iocc[plaq[0][p]] == 1 && lspin[plaq[4][p]] == 1) || (iocc[plaq[0][p]] == 0 && lspin[plaq[4][p]] == -1)){      
           continue;
        }
        if((iocc[plaq[0][p]] == 1 && lspin[plaq[4][p]] == -1) || (iocc[plaq[0][p]] == 0 && lspin[plaq[4][p]] == 1)){
          continue;
        } 
        else{ printf("error : invalid plaquette = %d, type2;\n", p); exit(1);}   
    
    }
  }


}



void fprintConf(FILE *fptr)
{
  int i, j, t;

  fprintf(fptr,"sites:\n");
  for(t=0; t<LT; t++)
  {
    for(j=0; j<LY; j++)
    {
      for(i=0; i<LX; i++)
      {
        fprintf(fptr,"%d ",iocc[i + LX * j + LX * LY * t]);
      }
    }
    fprintf(fptr,"\n");
  }
  fprintf(fptr,"links horizontal:\n");
  for(t=0; t<LT; t++)
  {
    for(j=0; j<LY; j++)
    {
      for(i=0; i<LX; i++)
      {
        fprintf(fptr,"%d ",lspin[i + LX * j + LX * LY * t]);
      }
    }
    fprintf(fptr,"\n");
  }
  fprintf(fptr,"links vertical:\n");
  for(t=0; t<LT; t++)
  {
    for(j=0; j<LY; j++)
    {
      for(i=0; i<LX; i++)
      {
        fprintf(fptr,"%d ",lspin[i + LX * j + LX * LY * t + LX*LY*LT]);
      }
    }
    fprintf(fptr,"\n");
  }
}
