//
//  zeros.h
//  hpc
//
//  Created by YIPENG ZHOU on 2018/3/14.
//  Copyright © 2018年 joe. All rights reserved.
//

#ifndef zeros_h
#define zeros_h

// zeros double array
void zerosda(int m, double* A){
    for (int i=0; i<m; i++){
        A[i] = 0.0;
    }
}

// zeros int array
void zerosia(int m, int* A){
    for (int i=0; i<m; i++){
        A[i] = 0;
    }
}

// zeros int matrix
void zerosim(int m, int n, int** A){
    for (int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            A[i][j] = 0;
        }
    }
}

// zeros double matrix
void zerosdm(int m, int n, double** A){
    for (int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            A[i][j] = 0.0;
        }
    }
}

#endif /* zeros_h */
