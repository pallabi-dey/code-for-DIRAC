#ifndef __DEFINE_H_INCLUDED__
#define __DEFINE_H_INCLUDED__

extern int LX,LY,LT,VOL,VOL2,V2,VOL4,V4,SVOL;
extern int HCLASS,SEED,ETERMS,NLOOPS,MERONNUM,MERONSECTORCOUNT,CONFIG,FLAG;
extern int histFlag;
/* fermion  and gauge variables */
extern int *iocc,*lspin;
/* variables to label lattice and neighbors */
extern int *ixc,*iyc,*itc;                 // serial to lattice coordinares
extern int *itup,*itdn, *ltup, *ltdn; // +t,-t neighbors for sites and links
extern int *xp1,*xm1;                      // +x neighbor
extern int *yp1,*ym1;                      // +y neighbor
extern int *parity;                        // parity
extern int *iin, *iout, *lout, *lin;  // spatial bonds in and out of a site
extern int *lsitel, *lsiter; //the sites that touch each link
extern int *ifwd,*ibwd;  //??not-used??
extern int *istag;
extern int *cflag,*lclus,*iclus;
extern int *breakup;
extern int STFLAG;
//int *lbc, *ltc;  //??not-used??
extern int arg, INDEX;
extern int itherm;
extern double beta,eps,t;
extern double prob, prob1,prob2;
// plaq, site and link defined in lattice.c navigate to fermions/links touching active plaquettes
extern int *plaq[6];
extern int *site[2];
extern int *link;
/*loop variables*/
extern int **loopLists, **loopLinkLists;
extern int *loopSizes, *loopLinkSizes;
extern int *whichLoop;
extern int loopNumber;
extern int alignment;

/* observables */
extern double eFlux, psiBpsi, fermN, eFlux0, psiBpsi0, plaqOp1, stagEfluxX, stagEfluxY;
extern double G2;
extern double *EEcorr, *CDW;
extern double EEcorrTot, CDWTot;
extern double *confHist;
extern int *sectorHist, *clusterHist, *iclusterHist, *lclusterHist;
extern int sector1, sector2, sector3;
/* structures for cluster operations */
typedef struct cluster
{
    int ID;
    int *p;
    int sizep;
    int meron;
    int *l;
    int sizel;
    struct cluster *prev;
    struct cluster *next;
} CLUSTER;
extern CLUSTER *firstcluster;
extern int totalclusters;
extern size_t idmax;

#endif

