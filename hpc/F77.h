//
//  F77.h
//  hpc
//
//  Created by YIPENG ZHOU on 2018/3/14.
//  Copyright © 2018年 joe. All rights reserved.
//

#ifndef F77_h
#define F77_h
#define F77NAME(x) x##_
extern "C"{
    // DNRM2 := sqrt( x'*x )
    void F77NAME(dnrm2)(const int& n, const double* x, const int& incx);
    //double a = F77NAME(dnrm2)(n,x,1);
    
    // 点积 a = vector x * vector y
    double F77NAME(ddot)(const int& n, const double* x, const int& incx,
                         const double* y, const int& incy);
    
    void F77NAME(dgemm)(const char& TransA,const char& TransB,
                        const int& m, const int& n,
                        const int& k, const double& alpha,
                        const double* A, const int& lda,
                        const double* B, const int& ldb,
                        const double& beta, double* C,
                        const int& ldc);
    
    void F77NAME(dgemv)(const char& TransA,
                        const int& m, const int& n,
                        const double& alpha,
                        const double* A, const int& lda,
                        const double* x, const int& incx,
                        const double& beta, double* y,
                        const int& incy);
    
    void F77NAME(dgbmv)(const char& TransA,
                        const int& m, const int& n,
                        const int& kl, const int& ku,
                        const double& alpha,
                        const double* A, const int& lda,
                        const double* x, const int& incx,
                        const double& beta, double* y,
                        const int& incy);
    
    void F77NAME(dgesv)(const int& N,const int& NRHS,double* A,
                        const int& LDA,int* IPIV,double* B,
                        const int& LDB, int& INFO );
    
    void F77NAME(dcopy)(const int& n, const double* x, const int& incx,
                        const double* y, const int& incy);
    
    
}

#endif /* F77_h */
