#include "Class_VTK_Wrappers.hpp"


void Demo_BasicVTK(){

    dvecarr3E    points ;
    ivector2D    connectivity ;
    
    dvector1D    pressure ;
    dvecarr3E    velocity ;
    
    
    points.resize(8) ;
    connectivity.resize(1) ;
    connectivity[0].resize(8) ;
    
    pressure.resize(8) ;
    velocity.resize(1) ;
    
    points[0][0]  =  0.  ;
    points[0][1]  =  0.  ;
    points[0][2]  =  0.  ;
    
    points[1][0]  =  1.  ;
    points[1][1]  =  0.  ;
    points[1][2]  =  0.  ;
    
    points[2][0]  =  0.  ;
    points[2][1]  =  1.  ;
    points[2][2]  =  0.  ;
    
    points[3][0]  =  1.  ;
    points[3][1]  =  1.  ;
    points[3][2]  =  0.  ;
    
    points[4][0]  =  0.  ;
    points[4][1]  =  0.  ;
    points[4][2]  =  1.  ;
    
    points[5][0]  =  1.  ;
    points[5][1]  =  0.  ;
    points[5][2]  =  1.  ;
    
    points[6][0]  =  0.  ;
    points[6][1]  =  1.  ;
    points[6][2]  =  1.  ;
    
    points[7][0]  =  1.  ;
    points[7][1]  =  1.  ;
    points[7][2]  =  1.  ;
    
    
    for( int i=0; i<8; i++)  connectivity[0][i] = i ; 
    
    for( int i=0; i<8; i++)  pressure[i] = (double) i ; 
    velocity[0][0] = 1. ;
    velocity[0][1] = 2. ;
    velocity[0][2] = 3. ;
    
    
//    { //Write only grid to VTK in ascii format
//        VtkUnstrVec   vtk(".", "ustr1", "ascii", 11, points, connectivity );
//        vtk.Write() ;
//    }
// 
//    { //Write grid and data to VTK in appended mode
//        VtkUnstrVec   vtk(".", "ustr2", "appended", 11, points, connectivity );
//        vtk.AddData( pressure, "press", "Point") ;
//        vtk.AddData( velocity, "vel", "Cell") ;
//        vtk.Write() ;
//    }

//    { //Read grid and data, rename and rexport
//        dvecarr3E    Ipoints ;
//        ivector2D    Iconnectivity ;
//        
//        dvector1D    Ipressure ;
//        dvecarr3E    Ivelocity ;
//    
// 
//        VtkUnstrVec   vtk(".", "ustr2", "appended", 11 );
//
//        vtk.AddData( Ipoints, "Points", "Point") ;
//        vtk.AddData( Iconnectivity, "connectivity", "Cell") ;
//        vtk.AddData( Ipressure, "press", "Point") ;
//        vtk.AddData( Ivelocity, "vel", "Cell") ;
//
//        vtk.Read() ;
//
//        Ipressure = 2. * Ipressure ;
//        
//        vtk.SetNames("./", "ustr3") ;
//        vtk.Write() ;
//    }

    { //Read grid and data from Paraview-generated file, rename and rexport
        vector< array<float,3>  >    Ipoints ;
        vector< vector<int64_t> >    Iconnectivity ;
        
        vector<float>   label ;
        vector<int64_t> ids ;
 
        VtkUnstrVec   vtk(".", "selection", "appended", 5 );

        vtk.AddData( Ipoints, "Points", "Point") ;
        vtk.AddData( Iconnectivity, "connectivity", "Cell") ;
        vtk.AddData( label, "STLSolidLabeling", "Cell") ;
        vtk.AddData( ids, "vtkOriginalCellIds", "Cell") ;

        vtk.Read() ;

        vtk.SetNames("./", "my_selection") ;
        vtk.Write() ;
    }
    

    return ;
};
