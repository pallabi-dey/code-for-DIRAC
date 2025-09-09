 #include<stdio.h>
 #include<math.h>
 #include<string.h>
 #include<stdlib.h>
 #include "ranlxd.h"
 #include "mkl.h"
 #include "define.h"

extern void printConf(void);
int charge(int);

/* Prints the gauge and fermion configurations together with the Gauss' Laws */
void printconf(void)
{
    int convert(int, int, int);
    int i, Q, j;
    for (i=0;i<VOL;i++) {
        j = xm1[i];
        Q = charge(i);
        printf("%d = (%d,%d,%d), (spin,ferm,spin) = (%d,%d,%d), G_x = %d\n", i, ixc[i], iyc[i], itc[i],
            lspin[j], iocc[i], lspin[i], Q);
    }
}

int charge(int site)
{
    int Q,exp; int len;
    int b1, b2, b3, b4;
    /*set star coords*/
    b1 = site;
    b2 = xm1[site];
    b3 = site+LX*LY*LT;
    b4 = ym1[site]+LX*LY*LT;
    
    /* Q_x = n_x - \sum_{i=1,2} [E_(x,x+e_i)-E_(x-e_i,x)] + [(-1)^{x_1+x_2}-1]/2  */
    exp = (ixc[site]+iyc[site])%2;

    if (exp==0) Q = 0;
    else   Q =-1;

    Q = Q + iocc[site] - (lspin[b1] - lspin[b2] + lspin[b3] - lspin[b4])/2;
    return Q;
}

void measure(void){
  extern int convert(int,int,int);
  int i;
  double b1, b2, b3, b4;
  double e1, e2, e3, e4;
  fermN = 0.0;
  eFlux = 0.0;
  psiBpsi = 0.0;
  stagEfluxX = 0.0;
  stagEfluxY = 0.0;
  for(i=0;i<VOL;i++){
    if(i < LX*LY)
      {
	eFlux0   = eFlux0 + .5 * lspin[i];
        if(LY > 1) eFlux0   = eFlux0 + .5 * lspin[i+LX*LY*LT];
        psiBpsi0 = psiBpsi0 + 2 * (.5 - istag[i])*iocc[i];
      }
    eFlux   = eFlux + .5 * lspin[i];
    if(LY > 1) eFlux   = eFlux + .5 * lspin[i+LX*LY*LT];
    psiBpsi = psiBpsi + (1 - 2*istag[i])*iocc[i];
    fermN   = fermN + iocc[i];
    stagEfluxX = stagEfluxX + 0.5*(1-2*istag[i])*lspin[i];
    stagEfluxY = stagEfluxY + 0.5*(1-2*istag[i])*lspin[i+LX*LY*LT];
  }
  eFlux   = eFlux/((double)LT);
  psiBpsi = psiBpsi/((double)LT);
  fermN   = fermN/((double)LT);
  stagEfluxX = stagEfluxX/((double)LT);
  stagEfluxY = stagEfluxY/((double)LT);
  
  plaqOp1 = 0.0;
  for(i=0;i<VOL;i++){
        e1 = 0.5*lspin[i];
        e2 = 0.5*lspin[xp1[i]+LX*LY*LT];
        e3 = 0.5*lspin[i+LX*LY*LT];
        e4 = 0.5*lspin[yp1[i]];
        if(e1>0)      plaqOp1 += fabs(((0.5+e1)*(0.5+e2)*(0.5-e3)*(0.5-e4)));
        else if(e1<0) plaqOp1 += fabs(((0.5-e1)*(0.5-e2)*(0.5+e3)*(0.5+e4)));
  }
  plaqOp1 /= ((double)LT); 
}


void configObs(FILE *fptr){
  int t,i,p;
  double e1,e2,e3,e4;
  double occ[SVOL];
  double xlink[SVOL];
  double ylink[SVOL];
  double plaqObs[SVOL];
  for(i=0;i<SVOL;i++){ occ[i] = 0.0; xlink[i] = 0.0; ylink[i] = 0.0; plaqObs[i]=0.0;}
  for(t=0;t<LT;t++){
      for(i=0;i<SVOL;i++){
          p = i+SVOL*t;
          e1 = 0.5*lspin[p];
          e2 = 0.5*lspin[xp1[p]+LX*LY*LT];
          e3 = 0.5*lspin[p+LX*LY*LT];
          e4 = 0.5*lspin[yp1[p]];
          occ[i]     += iocc[i+SVOL*t];
          xlink[i]   += 0.5*lspin[i+SVOL*t];
          ylink[i]   += 0.5*lspin[i+SVOL*t+LX*LY*LT];
          //plaqObs[i] += fabs(((0.5+e1)*(0.5+e2)*(0.5-e3)*(0.5-e4)));
          if(e1>0)      plaqObs[i] += fabs(((0.5+e1)*(0.5+e2)*(0.5-e3)*(0.5-e4)));
          else if(e1<0) plaqObs[i] += fabs(((0.5-e1)*(0.5-e2)*(0.5+e3)*(0.5+e4)));
          //printf("site=%d, xp1 %d yp1=%d, e1 = %le, e2 = %le, e3 = %le e4=%le\n", p, xp1[p], yp1[p], e1, e2, e3, e4);
      }     
  }
  for(i=0;i<SVOL;i++){
      occ[i]     /= ((double)LT);  
      xlink[i]   /= ((double)LT);
      ylink[i]   /= ((double)LT); 
      plaqObs[i] /= ((double)LT); 
      //printf("i = %d occ = %le xlink = %le ylink = %le flippable = %le\n", i, occ[i], xlink[i], ylink[i], plaqObs[i]);
      fprintf(fptr, "%d\t%.4le %.4le %.4le %.4le\n", i, occ[i], xlink[i], ylink[i], plaqObs[i]);
  }
}

void improvedMeasure(){
    int i,j;
    CLUSTER *tcluster;
    tcluster=firstcluster;
    psiBpsi = 0.0;
    eFlux = 0.0;
    for(i=0;i<totalclusters;i++){
        psiBpsi = 0.0;
        eFlux = 0.0;
        for(j=0;j<tcluster->sizep;j++){
            psiBpsi += istag[(tcluster->p)[j]]*iocc[(tcluster->p)[j]];
        }
        for(j=0;j<tcluster->sizel;j++){    
            eFlux   += 0.5*lspin[(tcluster->l)[j]];
        }
        tcluster=tcluster->next;  
        //printf("%le %le \n", psiBpsi, eFlux);   
    }
}

double GLsquared(int iter){
    double Gsq;
    int i;
    int Q, Qpls, Qmns;

    Gsq  = 0.0; 
    for(i=0;i<VOL;i++){
        Q = charge(i);
        Gsq += Q*Q;
    }
    Gsq = ((double)Gsq)/((double)LT);
    return Gsq;
}

void storeSector(void)
{
    int convert(int, int, int);
    int i, Q, j;
    sector1 = sector2 = sector3 = 0;
    int cond1 = 1; int cond2 = 1; 
    int cond3 = 1; int cond4 = 1;
    int cond5 = 1; int cond6 = 1;
    for (i=0;i<VOL;i++) {
        j = xm1[i];
        Q = charge(i);
        if(istag[i] == 0){ 
          if (Q != -1) cond1 = 0;
          if (Q !=  2) cond3 = 0;
          if (Q !=  1) cond5 = 0;
        }
        else{
          if (Q != 1) cond2 = 0;
          if (Q !=-2) cond4 = 0;
          if (Q !=-1) cond6 = 0;
        }  
      
          
    }
    if(cond1 && cond2) sector1 = 1;              // G = (-1,1) Sector
    if(cond3 && cond4) sector2 = 1;              // G = (2,-2) Sector
    if(cond5 && cond6) sector3 = 1;              // G = (1, 1) Sector
}



/*------------------------------------------------------------------------------------------------------------------------------------------------------------
double GLsquared(int iter) {
    double Gsq;
    int i;
    int Q, Qpls, Qmns;

    Gsq  = 0.0; Qpls = 0; Qmns = 0;
    for (i=0;i<VOL;i++) {
        Q = charge(i);
        Gsq += Q*Q;
        // collect the total positive charge contribution 
        if (Q >= 0) Qpls = Qpls + Q;
        else        Qmns = Qmns + Q;
        //printf("Qsquare is %d.\n",Q*Q);
    }
    // histogram the charge distribution 
    //if(Qpls>V2) printf("Qpls = %d\n",Qpls);
    if (Qpls < 0) {
        printf("Error. Qpls = %d\n", Qpls); exit(0);
    }
    if (Qmns > 0) {
        printf("Error. Qmns = %d\n", Qmns); exit(0);
    }
    if(ffabs(Qpls) != ffabs(Qmns)){
       //printf("Qpls = %d, Qmns = %d in iter = %d\n", Qpls, Qmns,iter);
       //printConf();
    }
    //Ghist[VOL+Qpls]++;
    //Ghist[VOL+Qmns]++;

    Gsq = ((double)Gsq)/((double)LT);
    return Gsq;
}

void allocSectorHist(void)
{
    int size,i;
    double gauss_law_dim;
    gauss_law_dim = 4.0;
    if(LY > 1) gauss_law_dim = 6.0;
    if(histFlag == 1)
    {
        size = (int) pow(gauss_law_dim,(double) LX*LY);
        printf("size is %d \n",size);
        sectorHist = mkl_calloc(size, sizeof(int), alignment);
    }
    else
    {
      sectorHist = mkl_calloc(2, sizeof(int), alignment);
    }

    //for(i=0; i<size; i++) printf("element %d of sectorHist is %d.\n",i,sectorHist[i]);
}

/*int updateMiniHist(void)
{
    size_t index, Q, ind1, ind2, b, ix, iy, i;
    int offset_e, offset_o, gauss_law_dim;

    offset_e = 1;
    offset_o = 2;
    gauss_law_dim = 4;
    if(LY > 1)
    {
        gauss_law_dim = 6;
        offset_e = 1;
        offset_o = 4;
    }

    index = 0;
    ind1 = 0;
    ind2 = 0;
    b = 1;
    for(ix=0; ix<LX; ix++)
    {
        for(iy=0; iy<LY; iy++)
        {
            if(((ix + iy)%2) == 0)
            {
                ind1 += offset_o * b;
                ind2 += offset_e * b;
                b *= gauss_law_dim;
            }
            else
            {
                ind1 += offset_e * b;
                ind2 += offset_o * b;
                b *= gauss_law_dim;
            }

        }
    }
    for(ix=0; ix<LX; ix++)
    {
        for(iy=0; iy<LY; iy++)
        {
            i=ix + LX * iy;
            Q = charge(i);
            //printf("site is %d, Q is %d\n", i, Q);
            if((ix+iy)%2 == 0) //even sites
            {
                // values are -1, 0, 1, 2
                index += (Q+offset_e) * (size_t) pow((double) gauss_law_dim, (double) i);
            }
            else //odd sites
            {
                // values are -2, -1, 0, 1
                index += (Q+offset_o) * (size_t) pow((double) gauss_law_dim, (double) i);
            }
        }
    }
    //printf("index is %zu\n",index);
    //printf("%zu %zu\n", ind1, ind2);

    if(index == ind2) sectorHist[0]++;
    else if(index == ind1) sectorHist[1]++;

    if(index == ind2)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int updateMiniHist(void)
{
    size_t Q, ix, iy, i;
    int type_1, type_2, max_flux;

    type_1 = 1;
    type_2 = 1;
    max_flux = 1;
    if(LY > 1) max_flux = 3;

  
    for(ix=0; ix<LX; ix++)
    {
        for(iy=0; iy<LY; iy++)
        {
            i=ix + LX * iy;
            Q = charge(i);
            //printf("site is %d, Q is %d\n", i, Q);
            if((ix+iy)%2 == 0) //even sites
            {
                if(Q != 0) type_1 = 0;
                if(Q != max_flux) type_2 = 0;
            }
            else //odd sites
            {
                if(Q != 0) type_1 = 0;
                if(Q != (-max_flux)) type_2 = 0;
            }
        }
    }
    //printf("index is %zu\n",index);
    //printf("%zu %zu\n", ind1, ind2);

    if(type_1 == 1) sectorHist[0]++;
    else if(type_2 == 1) sectorHist[1]++;

    if(type_1 == 1)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void updateSectorHist(void)
{
    int index, i, Q,parsite,offset_e,offset_o;
    double gauss_law_dim;
    offset_e = 1;
    offset_o = 2;
    gauss_law_dim = 4.0;
    if(LY > 1)
    {
        gauss_law_dim = 6.0;
        offset_e = 1;
        offset_o = 4;
    }
    // timeslice 0
    index = 0;
    for(i = 0; i<LX*LY; i++)
    {
        parsite = i;
        if(LY > 1) parsite = ixc[i]+iyc[i];
        Q = charge(i);
        //printf("site is %d, Q is %d\n", i, Q);
        if((parsite % 2) == 0) //even sites
        {
            // values are -1, 0, 1, 2, 3, 4
            //printf("Q is %d, %d\n", Q,2);
            index += (Q+offset_e) * (int) pow(gauss_law_dim, (double) i);
        }
        else //odd sites
        {
            // values are -4, -3, -2, -1, 0, 1
            //printf("Q is %d, %d\n", Q,4);
            index += (Q+offset_o) * (int) pow(gauss_law_dim, (double) i);
        }

    }
    //printf("index is %d\n",index);
    sectorHist[index] += 1;
}
*/
/*
void printSectorHist(FILE *fPtr)
{
    int i;
    int spa_vol;
    spa_vol = LX;
    if(LY > 1) spa_vol = LX*LY;
    if(histFlag == 1)
    {
        int size;
        size = (int) pow(6.0,(double) spa_vol);
        for(i=0; i < size; i++)
        {
            if(sectorHist[i] != 0)
            {
                fprintf(fPtr,"%d %d\n", i, sectorHist[i]);
            }*/
            /*else
            {
                printf("got zero\n");
            }*/
/*        }
    }
    else
    {
        for(i=0; i < 2; i++)
        {
            if(sectorHist[i] != 0)
            {
                fprintf(fPtr,"%d %d\n", i, sectorHist[i]);
            }
        }
    }

}

void printClustSizeHist(FILE *fptr)
{
    int i;
    for(i=0; i<2*VOL+1; i++)
    {
        if(clusterHist[i] != 0)
        {
            fprintf(fptr,"%d %d\n", i, clusterHist[i]);
        }
    }
}

void updateClustHist(void)
{
    CLUSTER *tcluster;
    tcluster = firstcluster;
    do
    {
        lclusterHist[tcluster->sizel]++;
        iclusterHist[tcluster->sizep]++;
        clusterHist[tcluster->sizel + tcluster->sizep]++;
        tcluster = tcluster->next;
    } while (tcluster != firstcluster);
}

/* The observables are:
 Eflux = <E> = (1/LT) \sum_{x,t} S^z (x,t)
 psiBarpsi   = (1/LT) \sum_{x,t} (-1)^x n (x,t)
 N_F         = (1/LT) \sum_{x,t} n (x,t)

 n_x = \psi^\dagger_x \psi_x

 EE-correlator = <E(r) E(0) > = (1/LT) \sum_{x,t} S^z (x,t) S^z (x+r,t)
 CDW-corr (r)  = (1/LT) \sum_{x,t} [ n(x,t) - 1/2 ] [ n(x+r,t) - 1/2 ] (-1)^r
 */

/* measure unimproved estimators */
/*
void measure(void){
  extern int convert(int,int,int);
  int i;
  //eFlux0 = 0.0; psiBpsi0 = 0.0;
  //eFlux = 0.0; psiBpsi = 0.0; 
  fermN = 0.0;
  for(i=0;i<VOL;i++){
    if(i < LX*LY)
      {
	    eFlux0   = eFlux0 + .5 * lspin[i];
        if(LY > 1) eFlux0   = eFlux0 + .5 * lspin[i+LX*LY*LT];
        psiBpsi0 = psiBpsi0 + 2 * (.5 - istag[i])*iocc[i];
      }
    eFlux   = eFlux + .5 * lspin[i];
    if(LY > 1) eFlux   = eFlux + .5 * lspin[i+LX*LY*LT];
    psiBpsi = psiBpsi + 2 * (.5 - istag[i])*iocc[i];
    fermN   = fermN + iocc[i];
  }
  eFlux   = eFlux/((double)LT);
  psiBpsi = psiBpsi/((double)LT);
  fermN   = fermN/((double)LT);
}
*/
/*improved estimators*/
/*
void improvedmeasure(void){
    CLUSTER *tcluster;
    int i,j,index;
    double cdwmer[2], eemer[2];
    tcluster=firstcluster;
    psiBpsi = 0.0;
    CDWTot = 0.0;
    eFlux = 0.0;
    EEcorrTot = 0.0;
    if(MERONSECTORCOUNT==0){
        for(j=0; j<totalclusters; j++)
        {
            //totalclusters
            psiBpsi = 0.0;
            eFlux = 0.0;
            for(i=0; i<tcluster->sizep; i++)
            {
                psiBpsi += istag[(tcluster->p)[i]]*(iocc[(tcluster->p)[i]]-0.5);
                eFlux += 0.5 * lspin[(tcluster->p)[i]];
            }
            CDWTot += psiBpsi * psiBpsi / ((double) VOL);
            EEcorrTot += eFlux * eFlux / ((double) VOL);
            tcluster=tcluster->next;
        }
    }
    else if(MERONSECTORCOUNT==2){
        cdwmer[0] = 0.0;
        cdwmer[1] = 0.0;
        eemer[0] = 0.0;
        eemer[1] = 0.0;
        index=0;
        for(j=0; j<totalclusters; j++)
        {
            //find merons
            if(tcluster->meron == 1)
            {
                for(i=0; i<tcluster->sizep; i++)
                {
                    cdwmer[index] += istag[(tcluster->p)[i]]*(iocc[(tcluster->p)[i]]-0.5);
                    eemer[index] += 0.5 * lspin[(tcluster->p)[i]];
                }
                index++;
                tcluster=tcluster->next;
            }
        }
        CDWTot += ffabs(cdwmer[0]) * ffabs(cdwmer[1]) * 2 / (pow(2.0,totalclusters)) / ((double) VOL);
        EEcorrTot += ffabs(eemer[0]) * ffabs(eemer[1]) * 2 / (pow(2.0,totalclusters)) / ((double) VOL);
    }
}
*/
/* measure correlation functions */
/*
void corrF(FILE *fptrC){
  extern int convert(int x,int y,int t);
  int x,t,r;
  int ind,indx;
  int sign;

  // initialize correlation functions 
  for(r=0; r<LX; r++){
     EEcorr[r]=0.0; CDW[r]=0.0;
  }
  EEcorrTot = 0.0, CDWTot = 0.0;

  // compute the correlation functions 
  sign=1;
  for(r=0; r<LX; r++)
  {
    for(t=0; t<LT; t++)
    {
        for(x=0; x<LX; x++)
        {
	        ind        = convert(x,0,t);
	        indx       = convert((x+r)%LX,0,t);
                EEcorr[r] += .25 * lspin[ind]*lspin[indx];
	        CDW[r]    += sign*(iocc[ind]-0.5)*(iocc[indx]-0.5);
        }
    }
    EEcorr[r] = EEcorr[r]/((double)LT);
    CDW[r]    = CDW[r]/((double)LT);
    sign = -sign;
  }

  for(r=0; r<LX; r++)
  {
    EEcorrTot += EEcorr[r];
    CDWTot += CDW[r];
    //fprintf(fptrC,"%d %lf %lf\n",r, EEcorr[r], CDW[r]);
  }
}*/
