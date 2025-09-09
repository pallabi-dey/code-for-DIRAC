#include <stddef.h>

int LX,LY,LT,VOL,VOL2,V2,VOL4,V4,SVOL;
int HCLASS,SEED,ETERMS,NLOOPS,MERONNUM,MERONSECTORCOUNT,CONFIG,FLAG;
int histFlag;
/* fermion  and gauge variables */
int *iocc,*lspin;
/* variables to label lattice and neighbors */
int *ixc,*iyc,*itc;                 // serial to lattice coordinares
int *itup,*itdn, *ltup, *ltdn; // +t,-t neighbors for sites and links
int *xp1,*xm1;                      // +x neighbor
int *yp1,*ym1;                      // +y neighbor
int *parity;                        // parity
int *iin, *iout, *lout, *lin;  // spatial bonds in and out of a site
int *lsitel, *lsiter; //the sites that touch each link
int *ifwd,*ibwd;  //??not-used??
int *istag;
int *cflag,*lclus,*iclus;
int *breakup;
int STFLAG;
//int *lbc, *ltc;  //??not-used??
int arg, INDEX;
int itherm;
double beta,eps,t;
double prob, prob1,prob2;
// plaq, site and link defined in lattice.c navigate to fermions/links touching active plaquettes
int *plaq[6];
int *site[2];
int *link;
/*loop variables*/
int **loopLists, **loopLinkLists;
int *loopSizes, *loopLinkSizes;
int *whichLoop;
int loopNumber;
int alignment;

/* observables */
double eFlux, psiBpsi, fermN, eFlux0, psiBpsi0, plaqOp1, stagEfluxX, stagEfluxY;
double G2;
double *EEcorr, *CDW;
double EEcorrTot, CDWTot;
double *confHist;
int *sectorHist, *clusterHist, *iclusterHist, *lclusterHist;
int sector1, sector2, sector3;
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
CLUSTER *firstcluster;
int totalclusters;
size_t idmax;

