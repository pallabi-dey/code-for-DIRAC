 #include<stdio.h>
 #include<math.h>
 #include<string.h>
 #include<stdlib.h>
 #include "ranlxd.h"
 #include "mkl.h"
 #include "define.h"
 
 
void readConf(){
    FILE *fptr = fopen("conf.bin", "rb");
    if(!fptr){
    perror("fopen");
    fprintf(stderr, "Error: Could not open equilibrium config.\n");
    exit(EXIT_FAILURE);
    }

    fread(iocc, sizeof(int), VOL, fptr);
    fread(lspin, sizeof(int), 2*VOL, fptr);
    fclose(fptr);
}


void writeConf() {
  //FILE *fptr = fopen("conf.bin", "wb");
  char filename[100];
  snprintf(filename, sizeof(filename), "conf-%d.bin", INDEX);
  FILE *fptr = fopen(filename, "wb");
  if (!fptr) {
    fprintf(stderr, "Error: Could not open file conf_equli.bin for writing\n");
    exit(EXIT_FAILURE);
  }

  if (fwrite(iocc, sizeof(int), VOL, fptr) != VOL) {
    fprintf(stderr, "Error: Failed to write iocc to conf_equli.bin\n");
    fclose(fptr);
    exit(EXIT_FAILURE);
  }

  if (fwrite(lspin, sizeof(int), 2 * VOL, fptr) != 2 * VOL) {
    fprintf(stderr, "Error: Failed to write lspin to conf_equli.bin\n");
    fclose(fptr);
    exit(EXIT_FAILURE);
  }

  fclose(fptr);
}


int writerand(const char *filename){
    FILE *fptr;
    int size = rlxd_size();
    int *state = (int *)malloc(size * sizeof(int));

    if(!state){
       fprintf(stderr, "Failed to allocate memory for RNG state\n");
       return -1;
    }

    rlxd_get(state);

    fptr = fopen(filename, "wb");
    if(!fptr){
        fprintf(stderr, "Failed to open file for writing: %s\n", filename);
        free(state);
        return -1;
    }

    if (fwrite(state, sizeof(int), size, fptr) != size) {
        fprintf(stderr, "Failed to write RNG state to file\n");
        fclose(fptr);
        free(state);
        return -1;
    }

    fclose(fptr);
    free(state);
    return 0;
}


int readrand(const char *filename){
    FILE *fptr;
    int size = rlxd_size();
    int *state = (int *)malloc(size * sizeof(int));

    if(!state){
        fprintf(stderr, "Failed to allocate memory for RNG state\n");
        return -1;
    }

    fptr = fopen(filename, "rb");
    if(!fptr){
        fprintf(stderr, "Failed to open file for reading: %s\n", filename);
        free(state);
        return -1;
    }

    if(fread(state, sizeof(int), size, fptr) != size){
        fprintf(stderr, "Failed to read RNG state from file\n");
        fclose(fptr);
        free(state);
        return -1;
    }

    fclose(fptr);
    rlxd_reset(state);
    free(state);
    return 0;
}
