#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include "ranlxd.h"
#include "mkl.h"
#include "define.h"

/* The plaq[0...5][0...VOL2-1] array keeps track of interactions encountered.
   The inner index runs as follows:

     3----5----2       0,1,2,3 denote the occupation numbers of fermions
     |         |
     |         |       4,5     denote the Sz state of the gauge links
     0----4----1
*/


int plaquettecategory(int p){
    // lspin, iocc
    /*
    o --- s --- o     e --- s --- e
    |           |     |           |
    |           |     |           |
    o --- s --- o     e --- s --- e
    */
    if(iocc[plaq[0][p]] == iocc[plaq[1][p]] && iocc[plaq[1][p]] == iocc[plaq[2][p]] && iocc[plaq[2][p]] == iocc[plaq[3][p]]){
       if(lspin[plaq[4][p]] != lspin[plaq[5][p]]){
          printf("error: invalid plaquette, type0!\n");
          exit(1);
       }
       return 0;
    }
     /*
         e --- d --- o     o --- u --- e
         |           |     |           |
         |           |     |           |
         o --- u --- e     e --- d --- o
     */
    else if(iocc[plaq[0][p]] != iocc[plaq[3][p]]){
       if(lspin[plaq[4][p]] != lspin[plaq[5][p]] && iocc[plaq[0][p]] == iocc[plaq[2][p]] && iocc[plaq[0][p]] != iocc[plaq[1][p]]){
          if((iocc[plaq[0][p]] == 1 && lspin[plaq[4][p]] == 1) || ((iocc[plaq[0][p]] == 0 && lspin[plaq[4][p]] == -1))){
          return 1; }
       }
       else{ printf("error : invalid plaquette, type1;\n");}   
     
    }
    
     /* A. (ONLY A Breakup)
         o --- d --- e     e --- u --- o
         |           |     |           |
         |           |     |           |
         o --- d --- e     e --- u --- o
         
       B. ( BOTH A && D Breakup)
         o --- u --- e     e --- d --- o
         |           |     |           |
         |           |     |           |
         o --- u --- e     e --- d --- o
     */
    else if(iocc[plaq[0][p]] == iocc[plaq[3][p]] && iocc[plaq[1][p]] == iocc[plaq[2][p]] && iocc[plaq[0][p]] != iocc[plaq[1][p]] 
            && lspin[plaq[4][p]] == lspin[plaq[5][p]]){
        if((iocc[plaq[0][p]] == 1 && lspin[plaq[4][p]] == -1) || (iocc[plaq[0][p]] == 0 && lspin[plaq[4][p]] == 1)){
           //printf("plaquette %d, (%d,%d,%d,%d), (%d,%d)\n",p, iocc[plaq[0][p]],iocc[plaq[1][p]],iocc[plaq[2][p]],
           //                        iocc[plaq[3][p]],lspin[plaq[4][p]],lspin[plaq[5][p]]);
           //printf("sites: %d %d %d %d %d %d\n",plaq[0][p],plaq[1][p],plaq[2][p],plaq[3][p],plaq[4][p],plaq[5][p]);       
           return 2;
        }
        if((iocc[plaq[0][p]] == 1 && lspin[plaq[4][p]] == 1) || (iocc[plaq[0][p]] == 0 && lspin[plaq[4][p]] == -1)){
           //printf("plaquette %d, (%d,%d,%d,%d), (%d,%d)\n",p, iocc[plaq[0][p]],iocc[plaq[1][p]],iocc[plaq[2][p]],
           //                         iocc[plaq[3][p]],lspin[plaq[4][p]],lspin[plaq[5][p]]);
           //printf("sites: %d %d %d %d %d %d\n",plaq[0][p],plaq[1][p],plaq[2][p],plaq[3][p],plaq[4][p],plaq[5][p]);
           return 3;
        } 
        else{ printf("error : invalid plaquette, type2;\n"); exit(1);}   
    
    }
    
    
}


/* l is the link, d is the direction to follow links (0 is up, 1 is down)*/
int next(int i, int *l, int *d){
   // breakup 0 is the A breakup
   //  |         |
   //  |         |
   // breakup 1 is the D breakup:
   //  -----------
   //
   //  -----------
   //printf("site argument is %d\n",i);
   //printf("breakup above is %d, breakup below is %d.\n",breakup[site[0][i]],breakup[site[1][i]]);
   //printf("cluster id's for sites are %d, %d, %d, %d.\n",iclus[itup[i]],iclus[iout[i]],iclus[itdn[i]],iclus[iin[i]]);
   
   if((breakup[site[0][i]] == 0 ) && iclus[itup[i]] == -1){ // if the breakup above i is 0 (A), we first see if we can move up
       *l = -1; //!added link value
       *d = -1; //!added link value
        return itup[i];
    }
    else if((breakup[site[0][i]] == 1) && iclus[iout[i]] == -1){ // if the breakup above i is 1 (D) or 2 (C),
                                                                 // we first see if we can go right (even) or left (odd)
        *l = plaq[4][site[0][i]];
        *d = 1;
        return iout[i];
    }
    else if((breakup[site[1][i]] == 0 ) && iclus[itdn[i]] == -1){ // if the breakup below i is 0 (A), 
                                                                  // and we can't go up or right (left), we then see if we can move down
        *l = -1;
        *d = -1;
        return itdn[i];
    }
    else if((breakup[site[1][i]] == 1 ) && iclus[iin[i]] == -1){ // if the breakup below i is 1 (D) or 2 (C),
                                                                 // and we can't go up or right (left), we then see if we can move left (right)
         *l = plaq[5][site[1][i]];
         *d = 0;
         return iin[i];
    }
    
    else{ // no next site--we've completed the loop
          *l = -1;
          *d = -1;
          return -1;
    }

}


void updateBreakup(){
   int i, j;
   double ran[1];
   int category;
   int breakup1[VOL2];
   int *list, *list1;
   int size, size1;
   int plaq_dir=0;
   
   list = NULL;
   size = 0;
   extern int plaquettecategory(int);
   extern void plaquetteClusters(int, int **, int *,int);
   extern int notin( int, int *, int );
   extern void Cluster(int *, int);
   
   /*printf("\n initial breakup list:\n");
   for(i=0;i<VOL2;i++){
       printf("%d ",breakup[i]);
   }
   printf("\n"); */
   
   for(i=0;i<VOL2;i++){
      category = plaquettecategory(i);
      if(category == 0) breakup1[i] = 0;
      else if(category == 1) breakup1[i] = 1;
      else if(category == 2) breakup1[i] = 0;
      else if(category == 3){
        if(ixc[plaq[0][i]] == ixc[plaq[1][i]]){ 
           prob = prob2;
           plaq_dir = 2; 
           //printf("plaqutte %d along y\n", i);
        }
        else if(iyc[plaq[0][i]] == iyc[plaq[1][i]]){ 
           prob = prob1;
           plaq_dir = 1;
           //printf("plaqutte %d along x\n", i);
        }  
        //if(plaq_dir == 1)      prob = prob1;
        //else if(plaq_dir == 2) prob = prob2;
        ranlxd(ran,1);
        if(ran[0] < prob) breakup1[i] = 0;
        else              breakup1[i] = 1;
      }
      
   }
   
   for(i=0;i<VOL2;i++){
       if(breakup1[i] != breakup[i]){
          breakup[i] = breakup1[i];
          list1 = NULL;
          size1 = 0;
          plaquetteClusters(i, &list1, &size1, breakup[i]);
          for(j=0;j<size1;j++){
             if(notin(list1[j], list, size) == 1) {
                size++;
                list = mkl_realloc(list, size * sizeof(int));
                list[size - 1] = list1[j];
             }
          }
          if(list1 != NULL){
             mkl_free(list1);
          }
          
       }   
   }
   /*printf("\n changed breakup list:\n");
   for(i=0;i<VOL2;i++){
       printf("%d ",breakup[i]);
   }
   printf("\n");
   printf("# of cluster will change =%d\n", size);
   printf("Cluster List are changed : ");
   for(int j=0;j<size;j++){
       printf("%d ", list[j]);
   }
   printf("\n");*/
   /* Start building cluster */
   Cluster(list, size);
   
   if (list != NULL) {
        mkl_free(list);
    }

}

int notin(int cID, int *cIDlist, int size){
   int i;
   if(cIDlist == NULL) return 1;
   else{
       for(i=0;i<size;i++){
         if(cID == cIDlist[i]) return 0;
       }
       return 1;
   }
}

/*find the clusters that touch a plaquette*/
void plaquetteClusters(int plaquette, int **clusterlist, int *nclusters, int newbreakup){
   int notin( int, int *, int );
   int i;
   int currlink;
   if(*clusterlist != NULL){
      mkl_free(*clusterlist);
      *clusterlist = NULL;
   }
   *nclusters = 0;
   //checking the sites
   for(i=0;i<4;i++){   
      if(notin(iclus[plaq[i][plaquette]], *clusterlist, *nclusters) == 1){
         (*nclusters)++;
         (*clusterlist) = mkl_realloc(*clusterlist, (*nclusters)*sizeof(int));
         (*clusterlist)[*nclusters-1] = iclus[plaq[i][plaquette]];
      }
   }
   //checking the links
   for(i=4;i<6;i++){  
      if(notin(lclus[plaq[i][plaquette]], *clusterlist, *nclusters) == 1){
         (*nclusters)++;
         (*clusterlist) = mkl_realloc(*clusterlist, (*nclusters)*sizeof(int));
         (*clusterlist)[*nclusters-1] = lclus[plaq[i][plaquette]];
      }
   }
}   

int matchID(CLUSTER *tcluster, int *clusterlist, int nclusters){
   int i;
   for(i=0;i<nclusters;i++){
      if(tcluster->ID == clusterlist[i]) return 1; //matches something on the list
   }
   return 0; //doesn't match anything
}



void Cluster(int *clusterlist, int nclusters){
   int i,j;
   extern int next(int, int *, int *);
   extern void addLink(int, int, int *, int **, int, int);
   extern void addSite(int, int, int *, int **, int, int);
   extern int matchID(CLUSTER *, int *, int);
   extern void allocateClusterSite(int, int , int **, int ***, int **, int ***);
   extern void allocateClusterLink(int, int , int **, int ***, int **, int ***);
   extern void deallocateCluster(int , int *, int **, int *, int **);
   extern void addCluster(CLUSTER *, CLUSTER *, int *, int, int *, int, int);
   extern void subtractCluster(CLUSTER *);
   extern void addLinkChain(int *, int, int, int *, int **, int, int);
   double ran[1];
   int currlink, currsite, direction;
   int *sitelist;
   int *linklist;
   int sites, links;
   int *clistsitesizes = NULL;
   int *clistlinksizes = NULL;
   int **clistsites = NULL;
   int **clistlinks = NULL;
   int *clistmerons = NULL;
   int clistoldmerons;
   int index=0;
   int newclusters = 0;
   int finishedloop;
   
   CLUSTER *tcluster;
   int moresites[VOL];
   int counterS = 0;
   int nextrasites = 0;
   sitelist = mkl_malloc(VOL * sizeof(int), alignment);
   linklist = mkl_malloc(2*VOL * sizeof(int), alignment);
   links = sites = 0;
   clistoldmerons = 0;
   int nW=0;
   int nH=0;
   int numberofloop;
   CLUSTER **uclusters;
   uclusters = mkl_malloc(nclusters * sizeof(CLUSTER *), alignment);
   for(i=0;i<nclusters;i++) uclusters[i] = NULL;
   for(i=0; i<VOL; i++) moresites[i] = -1;
   tcluster = firstcluster;
   /* Unmarked the sites of relevant clusters due to altered plaquette breakup */
   for(i=0;i<nclusters;i++){
      while(matchID(tcluster, clusterlist, nclusters) != 1){
            tcluster = tcluster->next;
      }
      uclusters[i] = tcluster;
      if(tcluster->meron == 1) clistoldmerons++;
      for(j=0;j<tcluster->sizep;j++){
         if((tcluster->p)[j] > VOL) printf("point is out of volume bounds: %d\n", (tcluster->p)[j]);
         sitelist[sites] = (tcluster->p)[j];
         sites++;
         iclus[(tcluster->p)[j]] = -1;
      }
      for(j=0;j<tcluster->sizel;j++){
         linklist[links] = (tcluster->l)[j];
         links++;
         lclus[(tcluster->l)[j]] = -1;
      }
      tcluster = tcluster->next;
   }
   /*printf("sites is %d\n",sites);
   for(i=0; i<sites; i++) printf("%d ",sitelist[i]);
   printf("\n");
   printf("links is %d\n",links);
   for(i=0; i<links; i++) printf("%d ",linklist[i]);
   printf("\n"); */
    
   /* Start building Clusters by searching through sitelist*/ 
   for(i=0;i<sites;i++){      
      if(iclus[sitelist[i]] == -1){   
         newclusters++;
         nextrasites = 0;
         counterS  = 0;
         nW = 0;
         nH = 0;
         currsite  = sitelist[i];
         allocateClusterSite(currsite, newclusters, &clistsitesizes, &clistsites, &clistlinksizes, &clistlinks);
         clistmerons = mkl_realloc(clistmerons, newclusters * sizeof(int));
         clistmerons[newclusters-1] = 0;  
         clistsites[newclusters - 1][0] = currsite;     
         iclus[currsite] = idmax+index; 
         finishedloop = 0;
         numberofloop = 1;
         //printf("starting site %d of cluster %ld\n", sitelist[i], idmax+index);
         
         /* 1. Add sites to a loop. 
            2. while connecting sites in the loop, include the links from D- breakup 
               and continue until reach another D-breeakup.
            3. Also, include the sites outside the current loop that are connected 
               to D-breakup (left & right sites)----> moresites.
            4. After completing the current loop and move to the moresites list 
               and add these sites to current loop and repeat from pt 1 until all 
               sites are exhausted.
         */       
         while(next(currsite, &currlink, &direction) != -1 || counterS<nextrasites){              
             currsite = next(currsite, &currlink, &direction);
             if(iclus[currsite] == -1){
                iclus[currsite] = idmax+index; 
                addSite(currsite, newclusters, clistsitesizes, clistsites, idmax, index); 
                //printf("site %d added to cluster %ld with nW = %d\n", currsite, idmax+index, nW);                 
                if(currlink != -1){
                  if(lclus[currlink] == -1){
                      addLink(currlink, newclusters, clistlinksizes, clistlinks, idmax, index);
                      nH++;
                      //printf("link %d added to cluster %ld\n", currlink, idmax+index);                      
                      if(direction == 0){
                         // cluster will go up and add the links until it encounters with D breakup plaq
                         currlink = ltup[currlink];
                         while(link[currlink] == -1){
                              addLink(currlink, newclusters, clistlinksizes, clistlinks, idmax, index);
                              //printf("Continued link %d added to cluster %ld\n", currlink, idmax+index); 
                              currlink = ltup[currlink];                              
                         }
                           
                       }
                       else if(direction == 1){
                         // cluster will go down and add the links until it encounters with D breakup plaq
                         currlink = ltdn[currlink];
                         while(link[currlink] == -1){
                              addLink(currlink, newclusters, clistlinksizes, clistlinks, idmax, index);
                              //printf("Continued link %d added to cluster %ld\n", currlink, idmax+index); 
                              currlink = ltdn[currlink];
                         }                      
                       }
                       /* for the sites connected to D-breakup */
                       if(breakup[link[currlink]] == 1 || breakup[link[currlink]] == 0){
                         addLinkChain(&currlink, direction, newclusters, clistlinksizes, clistlinks, idmax, index);
                         nH++;
                         if(iclus[lsitel[currlink]] == -1){ //site to the left of link
                            moresites[nextrasites] = lsitel[currlink];
                            nextrasites++;
                         }   
                         if(iclus[lsiter[currlink]] == -1){ // check +x neighbor
                            moresites[nextrasites] = lsiter[currlink];                         
                            nextrasites++;
                         }                                   
                       }                      
                  }                     
               }//Finish connecting the link to loop.
             }// Finish connecting the sites to the loop.
            if(next(currsite, &currlink, &direction) == -1){ 
               /* Go through the extra links (if there) with D-breakup */
               if(lclus[currlink]==-1 && direction != -1 ){
                  addLink(currlink, newclusters, clistlinksizes, clistlinks, idmax, index);
                  if(direction == 0){
                     currlink = ltup[currlink];
                     while(link[currlink] == -1){
                           addLink(currlink, newclusters, clistlinksizes, clistlinks, idmax, index);
                           //printf("added currlink after complete site loop = %d in cluster %ld\n", currlink, idmax+index);
                           currlink = ltup[currlink];
                     }
                  }
                  else if(direction == 1){
                       // cluster will go down and add the links until it encounters with D breakup plaq
                       //printf("currlink =%d in cluster %ld\n", currlink, idmax+index);
                       currlink = ltdn[currlink];
                       while(link[currlink] == -1){
                             addLink(currlink, newclusters, clistlinksizes, clistlinks, idmax, index);
                             //printf("added currlink =%d in cluster %ld\n", currlink, idmax+index);
                             currlink = ltdn[currlink];
                       }                     
                  }
               }  
               if(counterS<nextrasites){
                  //printf("more site %d added to cluster %ld\n", currsite, idmax+index); 
                  while(iclus[moresites[counterS]] != -1 && counterS < nextrasites -1){ //added extra sites
                       counterS++;
                  }   
                  currsite = moresites[counterS];
                  //if(iclus[currsite] == -1) addSite(currsite, newclusters, clistsitesizes, clistsites, idmax, index); 
                  
                  //printf("more currsite added = %d currlink = %d, direction=%d in cluster %ld\n", currsite, currlink, direction, idmax+index);    
                  if(iclus[moresites[counterS]] != -1) counterS++; //exhausted all the sites 
                  if(iclus[currsite] == -1){ numberofloop++;}
                    
              }     
            }     
         }//complete 1st while loop
         
         index++;      
         //printf("nW = %d, nH = %d, NLoop = %d counterS=%d, nextrasites=%d.\n", nW, nH, numberofloop, counterS, nextrasites);      
         
      } 
        
   }
   
   for(i=0;i<links;i++){
      if(lclus[linklist[i]] == -1){
         if(breakup[link[linklist[i]]] == 1){ printf("orphaned link found. link = %d\n", linklist[i]); exit(1);}
         else{
              currlink = linklist[i];
              newclusters++;
              allocateClusterLink(currsite, newclusters, &clistsitesizes, &clistsites, &clistlinksizes, &clistlinks);
              clistmerons = mkl_realloc(clistmerons, newclusters * sizeof(int));
              clistmerons[newclusters-1] = 0;  
              clistlinks[newclusters - 1][0] = currlink;
              lclus[linklist[i]] = idmax + index;
              //printf("with vertical breakup, orphaned link = %d with cluster %d \n",linklist[i], lclus[currlink]);
              for(j=0;j<LT-1;j++){
                    currlink = ltup[currlink]; 
                    //printf("currlink=%d\n",currlink);  
                    if(lclus[currlink] != -1){
                       printf("error: trying to assign same link to two different clusters. currlink =%d with %d\n", currlink, lclus[currlink]);
                       exit(1);
                    }
                    else{ 
                        addLink(currlink, newclusters, clistlinksizes, clistlinks, idmax, index);
                        //printf("added currlink =%d in cluster %ld\n", currlink, idmax+index); 
                    } 
                    
                 
              }
              index++;
              
         }
      }
   }
   /* Check if there are left any sites or links */
    for(i=0;i<sites;i++){ 
      if(iclus[sitelist[i]]==-1){
         printf(" sites still unmarked = %d\n", sitelist[i]); 
         exit(0);
      }
   }   
         
   for(i=0;i<links;i++){ 
       if(lclus[linklist[i]]==-1){
          printf(" links still unmarked = %d\n", linklist[i]); 
          exit(0);
       }
   }    
   
   /* Subtract the clusters of which altered plaquette breakup */
   for(i=0;i<nclusters;i++){
       subtractCluster(uclusters[i]);
   }
   //add new clusters--idmax is updated
   for(i=0;i<newclusters;i++){
       if(totalclusters == 0) addCluster(firstcluster, firstcluster, clistsites[i], clistsitesizes[i], clistlinks[i], clistlinksizes[i], clistmerons[i]);
       else addCluster(firstcluster->prev, firstcluster, clistsites[i], clistsitesizes[i], clistlinks[i], clistlinksizes[i], clistmerons[i]);
   }
   
   
   mkl_free(sitelist);
   mkl_free(linklist);
   deallocateCluster(newclusters, clistsitesizes, clistsites, clistlinksizes, clistlinks);
   mkl_free(uclusters);
}

void addLink(int currlink, int newclusters, int *clistlinksizes, int **clistlinks, int idmax, int index)
{
   clistlinksizes[newclusters-1]++;
   clistlinks[newclusters - 1] = mkl_realloc(clistlinks[newclusters - 1], clistlinksizes[newclusters-1] * sizeof(int) );
   lclus[currlink] = idmax + index;
   clistlinks[newclusters - 1][clistlinksizes[newclusters-1]-1] = currlink;
}

void addSite(int currsite, int newclusters, int *clistsitesizes, int **clistsites, int idmax, int index)
{
   clistsitesizes[newclusters-1]++;
   clistsites[newclusters - 1] = mkl_realloc(clistsites[newclusters-1], clistsitesizes[newclusters-1] * sizeof(int));
   iclus[currsite] = idmax + index;
   clistsites[newclusters - 1][clistsitesizes[newclusters-1] - 1] = currsite;
}

void updateConfig(void){
   void flipClusters(CLUSTER *);
   int i;
   double ran[1];
   CLUSTER *tcluster;
   tcluster = firstcluster;
   for(i=0;i<totalclusters;i++){
      ranlxd(ran,1);
      if(ran[0]<0.5) flipClusters(tcluster);
      tcluster = tcluster->next;
   }
}

void flipClusters(CLUSTER *tcluster){
   int i;
   //printf("sites:\n");
   for(i=0; i<tcluster->sizep; i++)
   {
      //printf("%d ", (tcluster->p)[i]); 
      iocc[(tcluster->p)[i]] = 1 - iocc[(tcluster->p)[i]];
   }
   //printf("\n");
   //printf("links:\n");
   for(i=0; i<tcluster->sizel; i++)
   {   
       //printf("%d ", (tcluster->l)[i]); 
       lspin[(tcluster->l)[i]] = -lspin[(tcluster->l)[i]];
   }
   //printf("\n");
   
}





void allocateClusterSite(int currsite, int newclusters, int **clistsitesizes, int ***clistsites, int **clistlinksizes, int ***clistlinks) {
   *clistsitesizes = mkl_realloc(*clistsitesizes, newclusters * sizeof(int));
   *clistlinksizes = mkl_realloc(*clistlinksizes, newclusters * sizeof(int));
   *clistsites = mkl_realloc(*clistsites, newclusters * sizeof(int *));
   *clistlinks = mkl_realloc(*clistlinks, newclusters * sizeof(int *));
   (*clistsitesizes)[newclusters-1] = 1;
   (*clistlinksizes)[newclusters-1] = 0;
   (*clistsites)[newclusters - 1] = NULL;
   (*clistlinks)[newclusters - 1] = NULL;
   (*clistsites)[newclusters - 1] = mkl_realloc((*clistsites)[newclusters-1], 1 * sizeof(int));
   //(*clistsites)[newclusters - 1][0] = currsite;
}

void allocateClusterLink(int currlink, int newclusters, int **clistsitesizes, int ***clistsites, int **clistlinksizes, int ***clistlinks) {
   *clistsitesizes = mkl_realloc(*clistsitesizes, newclusters * sizeof(int));
   *clistlinksizes = mkl_realloc(*clistlinksizes, newclusters * sizeof(int));
   *clistsites = mkl_realloc(*clistsites, newclusters * sizeof(int *));
   *clistlinks = mkl_realloc(*clistlinks, newclusters * sizeof(int *));
   (*clistsitesizes)[newclusters-1] = 0;
   (*clistlinksizes)[newclusters-1] = 1;
   (*clistsites)[newclusters - 1] = NULL;
   (*clistlinks)[newclusters - 1] = NULL;
   (*clistlinks)[newclusters - 1] = mkl_realloc((*clistlinks)[newclusters-1], 1 * sizeof(int));
   //(*clistlinks)[newclusters - 1][0] = currlink;
}

void deallocateCluster(int newclusters, int *clistsitesizes, int **clistsites, int *clistlinksizes, int **clistlinks) {
   int i;
   for(i = 0; i < newclusters; i++) {
      if(clistsites[i] != NULL) mkl_free(clistsites[i]);
      if(clistlinks[i] != NULL) mkl_free(clistlinks[i]);
   }
   mkl_free(clistsites);
   mkl_free(clistlinks);
   if(clistsitesizes != NULL) mkl_free(clistsitesizes);
   if(clistlinksizes != NULL) mkl_free(clistlinksizes);
}


void addLinkChain(int *currlink, int direction, int newclusters, int *clistlinksizes, int **clistlinks, int idmax, int index)
{
   int count=0;
   while(link[*currlink]==-1 || breakup[link[*currlink]] == 0)
   {
      count++;
      addLink(*currlink, newclusters, clistlinksizes, clistlinks, idmax, index);
      //may need to be more complicated in 2d
      if(direction == 0) *currlink = ltup[*currlink];
      else *currlink = ltdn[*currlink];
      if(count == LT+2)
      {
         printf("error: infinite loop.\n");
         exit(1);
      }
   }
   // add last link, coming from a second D-breakup or C-breakup
   addLink(*currlink, newclusters, clistlinksizes, clistlinks, idmax, index);
}

