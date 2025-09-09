#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include "ranlxd.h"
#include "mkl.h"
#include "define.h"


void addCluster(CLUSTER *icluster, CLUSTER *fcluster, int *sites, int nsites, int *links, int nlinks, int meron)
{
    CLUSTER *tcluster;
    int i;

    tcluster = (CLUSTER *)mkl_malloc(sizeof(CLUSTER), alignment);

    tcluster->meron = meron;
    tcluster->ID = idmax;
    if(nsites == 0)
    {
        tcluster->p = NULL;
    }
    else
    {
      tcluster->p = mkl_malloc(nsites * sizeof(int), alignment);
        for(i=0; i<nsites; i++)
        {
            (tcluster->p)[i] = sites[i];
            iclus[sites[i]] = idmax;
            if(sites[i] >= VOL)
            {
                printf("error: tried to add site %d which is not on the lattice.\n",sites[i]);
                exit(1);
            }
        }
    }
    tcluster->sizep = nsites;

    if(nlinks == 0)
    {
    tcluster->l = NULL;
    }
    else
    {
      tcluster->l = mkl_malloc(nlinks * sizeof(int), alignment);
        for(i=0; i<nlinks; i++)
        {
            (tcluster->l)[i] = links[i];
            lclus[links[i]] = idmax;
        }
    }
    tcluster->sizel = nlinks;


    if( icluster == NULL && fcluster == NULL)
    {
        firstcluster = tcluster;
        tcluster->prev = firstcluster;
        tcluster->next = firstcluster;
    }
    else
    {
        tcluster->prev = icluster;
        tcluster->next = fcluster;
        icluster->next = tcluster;
        fcluster->prev = tcluster;
    }
   idmax++;
   totalclusters++;
}

void subtractCluster(CLUSTER *scluster)
{
    CLUSTER *icluster, *fcluster;
    if(totalclusters == 0 || scluster == NULL)
    {
        printf("subbond error: no bonds to subtract.\n");
        exit(1);
    }

    if(totalclusters == 1)
    {
        firstcluster = NULL;
    }
    else
    {
        icluster = scluster->prev;
        fcluster = scluster->next;
        icluster->next = fcluster;
        fcluster->prev = icluster;
        if(scluster == firstcluster) firstcluster = icluster;
    }

    if(scluster->p != NULL) mkl_free(scluster->p);
    if(scluster->l != NULL) mkl_free(scluster->l);


    mkl_free(scluster);
    totalclusters--;
}

void freeClusterList(void)
{
    CLUSTER *tcluster;
    CLUSTER *icluster;
    CLUSTER *fcluster;
    if(firstcluster != NULL)
    {
        while(firstcluster->next != firstcluster)
        {
            subtractCluster(firstcluster);
        }

        if(firstcluster->p != NULL) mkl_free(firstcluster->p);
        if(firstcluster->l != NULL) mkl_free(firstcluster->l);
        mkl_free(firstcluster);

    }

}

void resetID(void)
{
    CLUSTER *tcluster;
    int i;
    int ID = 0;
    tcluster = firstcluster;
    if(tcluster != NULL)
    {
        do{
            tcluster->ID = ID;
            for(i=0; i<tcluster->sizep; i++) iclus[(tcluster->p)[i]] = ID;
            for(i=0; i<tcluster->sizel; i++) lclus[(tcluster->l)[i]] = ID;
            ID++;
            tcluster = tcluster->next;
        }while(tcluster != firstcluster);
    }
    idmax = ID;

    if( idmax != totalclusters)
    {
        printf("error: max id isn't number of clusters (resetID, clusterOps.c)\n");
        exit(1);
    }
}

void printClusters(void)
{
    CLUSTER *tcluster;
    tcluster = firstcluster;
    int i;
    printf("# of clusters = %d\n", totalclusters);
    do{
        printf("Cluster %d: %d sites and %d links.\n", tcluster->ID, tcluster->sizep, tcluster->sizel);
        if(tcluster->sizep > 0)
        {
            printf("sites: ");
            for(i=0; i<tcluster->sizep; i++)
                printf("%d ", (tcluster->p)[i]);
            printf("\n");
        }
        if (tcluster->sizel > 0)
        {
            printf("links: ");
            for (i=0; i<tcluster->sizel; i++)
                printf("%d ",(tcluster->l)[i]);
            printf("\n");
        }
        tcluster = tcluster->next;
    }while(tcluster != firstcluster);
}

void fprintClusters(FILE *fptr)
{
    CLUSTER *tcluster;
    tcluster = firstcluster;
    int i;
    do{
        fprintf(fptr,"Cluster %d: %d sites and %d links.\n", tcluster->ID, tcluster->sizep, tcluster->sizel);
        if(tcluster->sizep > 0)
        {
            fprintf(fptr,"sites: ");
            for(i=0; i<tcluster->sizep; i++)
                fprintf(fptr,"%d ", (tcluster->p)[i]);
            fprintf(fptr,"\n");
        }
        if (tcluster->sizel > 0)
        {
            fprintf(fptr,"links: ");
            for (i=0; i<tcluster->sizel; i++)
                fprintf(fptr,"%d ",(tcluster->l)[i]);
            fprintf(fptr,"\n");
        }
        tcluster = tcluster->next;
    }while(tcluster != firstcluster);
}

int computeTimeWindingNumber(CLUSTER *tcluster)
{
    int wind;
    int i;
    wind = 0;
    for(i=0; i<tcluster->sizep; i++)
    {
        if(itc[(tcluster->p)[i]] == LT-1 && (itc[(tcluster->p)[(i+1)%(tcluster->sizep)]] == 0))
            wind++;
        else if(itc[(tcluster->p)[i]] == 0 && (itc[(tcluster->p)[(i+1)%(tcluster->sizep)]] == LT-1))
            wind--;
    }
    return wind;
}

int computeXWindingNumber(CLUSTER *tcluster)
{
    int wind;
    int i;
    wind = 0;
    for(i=0; i<tcluster->sizep; i++)
    {
        if(ixc[(tcluster->p)[i]] == LX-1 && (ixc[(tcluster->p)[(i+1)%(tcluster->sizep)]] == 0))
            wind++;
        else if(ixc[(tcluster->p)[i]] == 0 && (ixc[(tcluster->p)[(i+1)%(tcluster->sizep)]] == LX-1))
            wind--;
    }
    return wind;
}

void fprintWindingNumbers(FILE *fptr)
{
    int computeTimeWindingNumber(CLUSTER *);
    int computeXWindingNumber(CLUSTER *);
    CLUSTER *tcluster;
    tcluster = firstcluster;
    int i;
    do{
        fprintf(fptr,"%d %d ",computeXWindingNumber(tcluster), computeTimeWindingNumber(tcluster));
        tcluster = tcluster->next;
    }while(tcluster != firstcluster);
    fprintf(fptr,"\n");
}
