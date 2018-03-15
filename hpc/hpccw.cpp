//
//  main.cpp
//  hpccw
//
//  Created by YIPENG ZHOU on 2018/3/5.
//  Copyright © 2018年 joe. All rights reserved.
//
#include <math.h>
#include <vector>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdlib.h>//atoi, strtol
#include <fstream>
#include "print.h"
#include "zeros.h"
#include "F77.h"
#include "gvtk.h"
using namespace std;



int main(int argc, const char * argv[]) {
    
    // defining materials conductivity
    const double kx = 250.0;
    const double ky = 250.0;
    const double kxy = 0.0;
    
    //define section
    //const double th = atof(argv[11]); //thickness [m]
    const double th = 0.2; //thickness [m]
    
//    // defining materials conductivity
//    const double kx = atof(argv[8]);
//    const double ky = atof(argv[9]);
//    const double kxy = atof(argv[10]);
//
//    //define section
//    const double th = atof(argv[11]); //thickness [m]
    
    //integration scheme
    const int gaussorder = 2;
    
    // defining geometry
    double a = 0.25;
    double h1 = 1.0; // height at left edge [m]
    double h2 = h1 * 1.3;  //height at beam right edge [m]
    double L = 3.0 * h1; //beam length [m]
    double b = -a * L + (h2 - h1) / L; //beam height constant
    int nelem_x = 15; //number of elements in x-direction
    int nelem_y = 8; //number of elements in y-direction
    
//    double a = atof(argv[2]);
//    double h1 = atof(argv[3]); // height at left edge [m]
//    double h2 = atof(argv[4]);  //height at beam right edge [m]
//    double L = atof(argv[5]); //beam length [m]
//    double b = -a * L + (h2 - h1) / L; //beam height constant
//    int nelem_x = atof(argv[6]); //number of elements in x-direction
//    int nelem_y = atof(argv[7]); //number of elements in y-direction
    
    
    int nnode_elem = 4; //number of nodes in each element
    int nNodeDof[4] = {1,1,1,1}; //number of Dof per node (1 = Temperature only)
    int total_Dofs = sizeof( nNodeDof ) / sizeof( nNodeDof[0] ) * nnode_elem;
    
    
    //calculation
    int nelem = nelem_x * nelem_y; // total number of elements
    int nnode = (nelem_x + 1) * (nelem_y + 1); // total number of nodes
    
    ///------------- calculation of Nodal coord matrix---------------------------
    
    double* x = new double[nelem_x+1];//alternative way: double x[nelem_x+1];
    for (int i = 0; i < nelem_x+1; i++){
        x[i] = i * L / nelem_x;
    }
    
    double* h = new double[nelem_x+1];
    for (int i = 0; i < nelem_x+1; i++){
        h[i] = a * pow(x[i],2) + b * x[i] + h1;
    }

    ///Y
    double** Y = new double*[(nelem_y+1)];
    for (int i=0; i<nelem_y+1; i++) Y[i]=new double[(nelem_x+1)];
    for (int rownr = 0; rownr < nelem_y+1; rownr++){
        for (int colnr = 0; colnr < nelem_x+1; colnr++){
            Y[rownr][colnr] = - h[colnr] / 2 + rownr * h[colnr] / nelem_y;
        }
    }
    
    //Coord
    double** Coord = new double*[nnode];
    for (int i=0; i<nnode; i++) Coord[i]=new double[2];
    for (int colnr = 0; colnr < nelem_x+1; colnr++){
        int k = 0;
        for (int i = (colnr) * (nelem_y + 1); i < (colnr + 1) * (nelem_y + 1); i++){
            Coord[i][0] = x[colnr];
            Coord[i][1] = Y[k][colnr];
            k++;
        }
    }

    ///NodeTopo
    int **NodeTopo = new int* [nelem_y+1];
    for (int i = 0; i<nelem_y+1; i++) NodeTopo[i] = new int[nelem_x+1];
    zerosim(nelem_y+1, 2, NodeTopo);
    
    int k =0;
    for (int colnr = 0; colnr < nelem_x+1; colnr++){
        for (int rownr = 0; rownr < nelem_y+1; rownr++){
            NodeTopo[rownr][colnr] = k;
            k++;
        }
    }
    print2dint(nelem_y+1, nelem_x+1, NodeTopo);
    
    
    
    ///------ElemNode, eNode, eCoord, ElemX, ElemY------
    int** ElemNode = new int*[nelem];
    for(int i=0; i<nelem; i++) ElemNode[i] = new int[5];
    
    int elemnr = 0;
    for (int colnr=0; colnr<nelem_x;colnr++){
        for (int rownr=0; rownr<nelem_y;rownr++){
            ElemNode[elemnr][0] = elemnr;
            ElemNode[elemnr][4] = NodeTopo[rownr + 1][colnr];  //# Lower left node
            ElemNode[elemnr][3] = NodeTopo[rownr + 1][colnr + 1];  //# Lower right node
            ElemNode[elemnr][2] = NodeTopo[rownr][colnr + 1];  //# Upper right node
            ElemNode[elemnr][1] = NodeTopo[rownr][colnr];  //# upper left node
            elemnr += 1;
        }
    }
    
    int* eNodes = new int[4];
    double* eCoord = new double[8];
    
    for(int i=0; i<nelem; ++i){
        // write eCoord in col major format
        int it=0;
        for(int k=0;k<2;k++){
            for(int j=0; j<4; j++){
                eNodes[j] = ElemNode[i][j+1];
                eCoord[it] = Coord[eNodes[j]][k];
                it++;
            }
        }
    }
    
    
    
    ///------global dof-------
    int** globDof = new int*[nnode];
    for(int i=0; i<nnode; i++) globDof[i] = new int[2];
    zerosim(nnode, 2, globDof);
    
    int* globNodes = new int[nelem];
    int nNode = 0;
    int nDof = 0;
    int eDof = 0;
    
    for(int j=0; j<nelem; j++){
        
        // Global node numbers of element nodess
        for(int a=0;a<4;a++){
            globNodes[j] = ElemNode[j][a+1];
        }
        
        //loop over element nodes creates the first column
        for(int k=0;k<nnode_elem;k++){
            nNode=ElemNode[j][k+1];
            
            /* if the already existing ndof of the present node is less than
             the present elements ndof then replace the ndof for that node*/
            if(globDof[nNode][0]<nNodeDof[k]){
                globDof[nNode][0]=nNodeDof[k];
            }
            
        }
        
    }

    // counting the global dofs and inserting in globDof
    for(int j=0;j<nnode;j++){
        eDof = globDof[j][0];
        
        for(int k=0;k<eDof;k++){
            globDof[j][k+1] = nDof;
            nDof += 1;
        }
    }
    
    
    
///------ Assembly of global stiffness matrix K ------
    int gauss = gaussorder; // Gauss order
    double GP[2] = {-1 / sqrt(3),1 / sqrt(3)}; // Points
    
    // Weights
    int* W = new int[2];
    W[0] = 1;
    W[1] = 1;
    
    
    //conductivity matrix D
    double* D = new double[4];
    D[0] = kx;
    D[1] = D[2] = kxy;
    D[3] = ky;
    
    double** K = new double*[nDof];
    for(int i=0; i<nDof; i++) K[i] = new double[nDof];
    zerosdm(nDof, nDof, K);
    
    for(int i=0; i<nelem; ++i){
        // write eCoord in col major format
        int it=0;
        for(int k=0;k<2;k++){
            for(int j=0; j<4; j++){
                eNodes[j] = ElemNode[i][j+1];
                eCoord[it] = Coord[eNodes[j]][k];
                //cout<<eCoord[it] <<endl;
                it++;
            }
        }
        
        int* gDof = new int[nnode_elem];
        zerosia(nnode_elem, gDof);
        
        //produce gDof
        for(int j=0; j<nnode_elem;j++){
            for(int k=1;k < (nNodeDof[j]+1); k++){
                gDof[j] = globDof[eNodes[j]][k];
            }
        }
        
        
        //local stiffness matrix Ke by Gauss integration
        double* Ke = new double[nnode_elem*nnode_elem];
        zerosda(nnode_elem*nnode_elem, Ke);
        double* J = new double[8];
        double* DetJ = new double[gauss*gauss]; //for storing the dterminants fo J
        double* BDB = new double[16];
        double* N = new double[4];
        double* GN = new double[8];
        zerosda(8, GN);
        double DJ = 0.0;
        
        int info = 0;
        int* ipiv = new int[2];
        
        for(int i=0; i<gauss; i++){
            for(int j=0; j<gauss; j++){
                double eta = GP[i];
                double xi = GP[j];
                
                // SHAPE FUNCTIONS MATRIX
                N[0]=(1-xi)*(1-eta)/4;
                N[1]=(1+xi)*(1+eta)/4;
                N[2]=(1+xi)*(1-eta)/4;
                N[3]=(1-xi)*(1+eta)/4;
                
                // derivative (Gradient) of the shape functions
                GN[0]=-(1-eta)/4;
                GN[1]=(-(1-xi))/4;
                GN[2]=(1-eta)/4;
                GN[3]=(-(1+xi))/4;
                GN[4]=(1+eta)/4;
                GN[5]=(1+xi)/4;
                GN[6]=(-(1+eta))/4;
                GN[7]= (1-xi)/4;
                
                F77NAME(dgemm)('N','N',2,2,4,1.0,GN,2,eCoord,4,0.0,J,2); // Jacobian Matrix
                //printda(8,J);
                DJ = J[0]*J[3] - J[2]*J[1];
                //cout<<DJ<<endl;
                F77NAME(dgesv)(2,4, J, 2, ipiv, GN, 2, info); // matrix B is output GN
                //printMatrix(2, 4, GN, 2);
                F77NAME(dgemm)('T','N',4,2,2,1.0,GN,2,D,2,0.0,J,4); // output J is B_t*D
                //printMatrix(4, 2, J, 4);
                F77NAME(dgemm)('N','N',4,4,2,1.0,J,4,GN,2,0.0,BDB,4); // BDB is B_t * D * B
                //printMatrix(4, 4, BDB, 5);
                
                for (int k=0; k<4*4; k++){
                    Ke[k] = Ke[k] + BDB[k] * th * DJ * W[i] * W[j];
                }
            }
        }
        
        // input 1D array of ke into 2D array form KE
        double** KE = new double*[4];
        for(int i=0; i<4; i++) KE[i] = new double[4];
        
        for(int i=0; i<4; i++){
            for(int j=0; j<4; j++){
                KE[i][j] = Ke[j*4+i];
            }
        }
        
        // "inserting" the global element stiffness matrix into the global system stiffness matrix
        for (int i=0; i<4; i++){
            for(int j=0; j<4; j++){
                K[gDof[i]][gDof[j]]=K[gDof[i]][gDof[j]]+KE[i][j];
            }
        }
    }
    
    
    /// ------ SWITCHING CASES ------
    int T0;
    int q;
    int nTempNodes;
    int nFluxNodes;
    int lenTE;
    int* TempNodes = new int[nDof];
    int* fluxNodes = new int[nDof];
    int casenumber = atoi(argv[1]);
    
    switch (casenumber) {
        case 1:
            
            // Temperature B.C. left edge
            zerosia(nelem_y+1, TempNodes);
            nTempNodes = nelem_y+1; // Number of nodes with temp BC
            lenTE = nTempNodes;
            T0 = 10; // Temperature at boundary
            for(int i=0; i<nelem_y+1; i++){TempNodes[i] = NodeTopo[i][0];} // Nodes at the left edge of the beam
            
            // Define edges
            q = 2500; // constant flux at right edge of the beam
            nFluxNodes = nelem_y+1;
            for(int i=0; i<(nelem_y+1); i++) {
                fluxNodes[i] = NodeTopo[i][nelem_x];
            }
            
            break;
            
        case 2:
            
            // Temperature B.C. left edge
            zerosia(nelem_x+1, TempNodes);
            T0 = 10; // Temperature at boundary
            nTempNodes = nelem_x+1; // Number of nodes with temp BC
            lenTE = nTempNodes;
            for(int i=0; i<nelem_x+1; i++){TempNodes[i] = NodeTopo[0][i];} // Nodes at the left edge of the beam
            
            // Define edges
            q = 2500; // constant flux at right edge of the beam
            nFluxNodes = nelem_x+1;
            for(int i=0; i<(nelem_x+1); i++) {fluxNodes[i] = NodeTopo[nelem_y][i];}
            
            break;
            
        case 3:
            
            // Temperature B.C. left edge
            zerosia(nelem_y+1, TempNodes);
            T0 = -20; // Temperature at boundary
            nTempNodes = nelem_y+1; // Number of nodes with temp BC
            lenTE = nTempNodes;
            for(int i=0; i<nelem_y+1; i++){TempNodes[i] = NodeTopo[i][0];} // Nodes at the left edge of the beam
            
            // Define edges
            q = -5000; // constant flux at right edge of the beam
            nFluxNodes = nelem_x+1;
            for(int i=0; i<(nelem_x+1); i++) {fluxNodes[i] = NodeTopo[0][i];}
            
            
            break;
            
        default:
            cout << "retype case number!" << endl;
            break;
    }
    
    
    
    /// ------- after CASEs -------
    
    int** n_bc = new int*[4];
    for(int i=0; i<4; i++) n_bc[i] = new int[nFluxNodes-1];
    zerosim(4, (nFluxNodes-1), n_bc);
    
    int* ones = new int[nFluxNodes-1];
    for(int i=0; i<nFluxNodes-1; i++){ones[i] = 1;}
    
    for(int i=0; i<(nFluxNodes-1); i++){
        n_bc[0][i] = fluxNodes[i];
        n_bc[2][i] = q * ones[i]; // flux value at node 1
        n_bc[3][i] = q * ones[i]; // flux value at node 2
    }
    for(int i=1; i<nFluxNodes; i++){
        n_bc[1][i-1] = fluxNodes[i];
    }
    
    int nbe = nFluxNodes - 1; // Number of elements with flux load
    
    // Nodal Coordinates
    double* Coordx = new double[nnode];
    double* Coordy = new double[nnode];
    for(int i=0; i<nnode; i++){
        Coordx[i] = Coord[i][0];
        Coordy[i] = Coord[i][1];
    }
    
    double* f = new double[nDof]; // initialize nodal flux vector
    int diff = nFluxNodes - 4;
    int* n_bce = new int[diff];
    
    for(int i=0; i < nbe; i++){
        double* fq = new double[2];
        zerosda(2, fq); // initialize the nodal source vector
        int node1 = n_bc[0][i];
        int node2 = n_bc[1][i];
        
        for(int k=2; k<4; k++){
            n_bce[k-2] = n_bc[k][i];
        }
        
        double x1 = Coordx[node1]; // x coord of the first node
        double y1 = Coordy[node1]; // y coord of the first node
        double x2 = Coordx[node2]; // x coord of the second node
        double y2 = Coordy[node2]; // y coord of the second node
        
        double leng = sqrt(pow((x2-x1),2) + pow((y2-y1),2)); // edge length
        double detJ = leng / 2; // 1D Jacobian
        
        // integrate in xi direction (1D integration)
        for(int i=0; i<gauss; i++){
            double xi = GP[i];
            // 1D shape functions in parent domain
            double* N = new double[2];
            N[0] = 0.5 * (1-xi);
            N[1] = 0.5 * (1+xi);
            double flux = N[0]*n_bce[1]+N[1]*n_bce[0];
            for(int k=0; k<2; k++){
                fq[k] = fq[k] + W[i] * N[k] * flux * detJ * th; // nodal flux
            }
        }
        
        // define flux as negative integrals
        for(int k=0; k<2; k++){
            fq[k] = -fq[k];
        }
        
        f[node1] += fq[0];
        f[node2] += fq[1];
    }
    //printda(nDof, f);
    /// ------ Apply boundary conditions -------
    
    int** BC = new int*[nTempNodes];
    for(int i=0; i<nTempNodes; i++){BC[i]=new int[2];}
    
    for (int i=0; i<nTempNodes; i++){
        BC[i][0] = TempNodes[i];
        BC[i][1] = T0;
    }
    //print2dint(nTempNodes, 2, BC);
    
    /// ------- Assembling global "Force" vector -------
    int* OrgDof = new int[nDof];
    zerosia(nDof, OrgDof);
    
    double* T = new double[nDof]; //initialize nodal temperature vector
    zerosda(nDof, T);
    int rDof = nDof; // Reduced number of DOF
    
    int* ind = new int[nTempNodes];
    for(int i=0; i<nTempNodes; i++){ind[i] = BC[i][0];}
    
    for(int i=0; i<nTempNodes; i++){
        OrgDof[ind[i]] = -1;
    }
    
    for(int i=0; i<nTempNodes; i++){
        T[ind[i]] = BC[i][1];
    }
    
    rDof = rDof - nTempNodes;
 
    int* RedDof = new int[rDof];
    int counter1 = 0;
    for(int j=0; j<nDof; j++){
        if (OrgDof[j] == 0)
            OrgDof[j] = counter1;
        RedDof[counter1] = j;
        counter1 += 1;
    }
    
    /// ------ Partition Matrix ------
    // Init the array here
    int* mask_E = new int[nnode];
    zerosia(nnode, mask_E);
    int* invmask_E = new int[nnode];
    double* T_E = new double[lenTE];
    zerosda(lenTE, T_E);
    
    
    // mask_E
    int itT_E = 0;
    mask_E[0] = 1;
    
    
    for(int i=1; i<nDof; i++){
        if(T[TempNodes[i]]!=0){
            mask_E[TempNodes[i]] = 1;
            T_E[itT_E] = T[TempNodes[i]];
            itT_E ++;
        }
    }
    
    // ~mask_E
    int lenfF = nDof - lenTE;
    int it3 = 0;
    int itf_F = 0;
    double* f_F = new double[lenfF];
    for(int i=0; i<nnode; i++){
        if(mask_E[i]==1){
            invmask_E[i] = 0;
        }
        else if (mask_E[i]==0){
            invmask_E[it3] = i;
            f_F[itf_F] = f[i];
            itf_F ++;
            it3 ++;
        }
    }
    
    // define K_EE
    double* K_EE = new double[lenTE * lenTE];
    int itKEE = 0;
    for(int i=0; i<lenTE; i++){
        for(int j=0; j<lenTE; j++){
            K_EE[itKEE] = K[TempNodes[i]][TempNodes[j]];
            itKEE++;
        }
    }
    

    // define K_FF
    double* K_FF = new double[lenfF * lenfF];
    zerosda(lenfF * lenfF, K_FF);
    
    int itKFF = 0;
    for(int j=0; j<lenfF; j++){
        for(int i=0; i<lenfF; i++){
            K_FF[itKFF] = K[invmask_E[i]][invmask_E[j]];
            itKFF++;
        }
    }
    
    // define K_EF (store in ColMajor)
    double* K_EF = new double[lenTE * lenfF];
    zerosda(lenTE * lenfF, K_EF);
    
    int itKEF = 0;
    for(int j=0; j<lenfF; j++){
        for(int i=0; i<lenTE; i++){
            K_EF[itKEF] = K[TempNodes[i]][invmask_E[j]];
            itKEF++;
        }
    }
    
    /// ------ solve for d_F ------
    double* y = new double[lenfF];
    zerosda(lenfF, y);
    F77NAME(dgemv)('T',lenTE,lenfF,1.0,K_EF,lenTE,T_E,1,0.0,y,1);
    
    double* rhs = new double[lenfF];
    zerosda(lenfF, rhs);
    for(int i=0; i<lenfF; i++){
        rhs[i] = f_F[i] - y[i];
    }
    
    int info1 = 0;
    int* ipiv1 = new int[lenfF];
    F77NAME(dgesv)(lenfF, 1, K_FF, lenfF, ipiv1, rhs, lenfF, info1);
    
    double* T_F = new double[lenfF];
    zerosda(lenfF, T_F);
    
    F77NAME(dcopy)(lenfF, rhs, 1, T_F, 1); // copy solution to T_F
    
    /// ------ reconstruct the global displacement d -------
    for(int i=0; i<lenTE; i++){
        T[TempNodes[i]] = T_E[i];
    }
   
    for(int i=0; i<lenfF; i++){
        T[invmask_E[i]] = T_F[i];
    }
  
    /// ------ compute the reaction f_E ------
    double* f_E = new double[lenTE];
    zerosda(lenTE, f_E);
    
    double* s1 = new double[lenTE];
    zerosda(lenTE, s1);
    double* s2 = new double[lenTE];
    zerosda(lenTE, s2);
    
    F77NAME(dgemv)('N', lenTE, lenTE, 1.0, K_EE, lenTE, T_E, 1, 0.0, s1, 1);
    F77NAME(dgemv)('N', lenTE, lenfF, 1.0, K_EF, lenTE, T_F, 1, 0.0, s2, 1);
    
    for(int i=0; i<lenTE; i++){
        f_E[i] = s1[i] + s2[i];
    }
    
    // reconstruct the global reactions f
    for(int i=0; i<lenTE; i++){
        f[TempNodes[i]] = f_E[i];
    }
    for(int i=0; i<lenfF; i++){
        f[invmask_E[i]] = f_F[i];
    }
    printda(nDof, f);
    //generateVtk(Coord, nnode, nelem,nelem_x,nelem_y, ElemNode, T);
    
}











