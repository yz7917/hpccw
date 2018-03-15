//
//  gvtk.h
//  hpc
//
//  Created by YIPENG ZHOU on 2018/3/15.
//  Copyright © 2018年 joe. All rights reserved.
//

#ifndef gvtk_h
#define gvtk_h
/// generate vtk file
void generateVtk(double** Coord, int nnode, int nelem, int nelem_x, int nelem_y, int** ElemNode, double* T){
    
    ofstream vtkfile; //创建个对象
    vtkfile.open("vtk.txt");
    
    vtkfile <<"vtk output" << endl;
    vtkfile <<"ASCII" << endl;
    vtkfile <<"DATASET UNSTRUCTURED_GRID" << endl;
    vtkfile <<"POINTS 66 DOUBLE" << endl;
    
    // Coord
    for(int i=0;i<nnode;i++){
        for(int j=0;j<2;j++){
            vtkfile << Coord[i][j] <<" ";
        }
        vtkfile << "0.0" <<" ";
    }
    cout << endl;
    
    vtkfile <<"CELLS" << " " << nelem <<" "<< nelem * nelem << endl;
    
    // ElemNode
    for(int i=0;i<nelem;i++){
        vtkfile << ElemNode[i][0] << " " << ElemNode[i][1] << " " << ElemNode[i][2] << " " << ElemNode[i][3] << " " << ElemNode[i][4] << endl;
    }
    
    //
    vtkfile <<"CELL_TYPES" << nelem << endl;
    for(int i=0;i<nelem;i++){
        vtkfile << "9" << endl;
    }
    
    vtkfile <<"POINT_DATA" << nnode << endl;
    vtkfile <<"FIELD FieldData 1"<< endl;
    vtkfile <<"disp"<< " "<< nnode <<" "<< "double" << endl;
    for(int i=0; i<nnode; i++){
        vtkfile << T[i] <<" ";
    }
    
    vtkfile.close(); //关闭
    
    
}

#endif /* gvtk_h */
