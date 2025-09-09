#include <stdio.h>
#include <stdlib.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <bitset>
#include <vector>

typedef Eigen::SparseMatrix<double> SparseMatrix;
void writeH(char *fname, Eigen::SparseMatrix<double> &Hmatrix){
    int row, col, nonZero;
    int r, c, k;
    double val;
    
    FILE *fptr = fopen(fname, "w");
    if(!fptr){
      printf("could not open file %s to write\n", fname);
      exit(1);
    }

    row     = Hmatrix.rows();
    col     = Hmatrix.cols();
    nonZero = Hmatrix.nonZeros();

    fwrite(&row, sizeof(int), 1, fptr);
    fwrite(&col, sizeof(int), 1, fptr);
    fwrite(&nonZero, sizeof(int), 1, fptr);

    for(k=0;k<Hmatrix.outerSize();k++){
        for(SparseMatrix::InnerIterator it(Hmatrix, k); it; ++it){
            r   = it.row();
            c   = it.col();
            val = it.value();

            fwrite(&r, sizeof(int), 1, fptr);
            fwrite(&c, sizeof(int), 1, fptr);
            fwrite(&val, sizeof(double), 1, fptr);
        }
    }

    fclose(fptr);
    printf("Sparse matrix saved to %s\n", fname);
}


void readH(char *fname, Eigen::SparseMatrix<double> &Hmatrix){
    int row, col, nonZero;
    int r, c, k;
    double val;
    
    FILE *fptr = fopen(fname, "r");
    if(!fptr){
      printf("could not open file %s to write\n", fname);
      exit(1);
    }

    
    fread(&row, sizeof(int), 1, fptr);
    fread(&col, sizeof(int), 1, fptr);
    fread(&nonZero, sizeof(int), 1, fptr);

    Hmatrix.resize(row, col);

    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nonZero);

 
    for(k=0;k<nonZero;k++){
        fread(&row, sizeof(int), 1, fptr);
        fread(&col, sizeof(int), 1, fptr);
        fread(&val, sizeof(double), 1, fptr);

        tripletList.emplace_back(row, col, val);
    }

    Hmatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    fclose(fptr);
    printf("read the Sparse matrix.\n");
}

