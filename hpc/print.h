//
//  print.h
//  hpc
//
//  Created by YIPENG ZHOU on 2018/3/14.
//  Copyright © 2018年 joe. All rights reserved.
//
#include<iostream>
#include <iomanip>
using namespace std;
#ifndef print_h
#define print_h

/// output 2D arrary double
void print2ddouble(int m, int n, double **Y){
    for (int rownr = 0; rownr < m; rownr++){
        for (int colnr = 0; colnr < n; colnr++){
            cout << setprecision(3)<<Y[rownr][colnr]<< " ";
        }
        cout << endl;
    }
}

// output 2D array int
void print2dint(int m, int n, int **Y){
    for (int rownr = 0; rownr < m; rownr++){
        for (int colnr = 0; colnr < n; colnr++){
            cout << setprecision(3)<<Y[rownr][colnr]<< " ";
        }
        cout << endl;
    }
}

// output int array
void printia(int m, int *Y){
    for (int rownr = 0; rownr < m; rownr++){
        cout << setprecision(3)<<Y[rownr] << endl;
    }
}

// output double array
void printda(int m, double *Y){
    for (int rownr = 0; rownr < m; rownr++){
        cout << setprecision(3)<<Y[rownr] << endl;
    }
}

// output colMajor vector
void printVector(int n, double* v, int prec){
    for (int i = 0; i < n; i++){
        cout  << setprecision(prec) << setw(5)<< v[i] << endl;
    }
}

///] output colMajor matrix
void printMatrix(int m, int n, double* Y, int prec){
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            cout  << setprecision(prec) << Y[j*m+i] << " ";
        }
        cout << endl;
    }
}


#endif /* print_h */
