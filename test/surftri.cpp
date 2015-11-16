// ========================================================================== //
//                - Class_SurfTri - Examples of usage                         //
//                                                                            //
// Example of usage for Class_SurfTri (grid manager for unstructured meshes)  //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author   : Alessandro Alaia                                                //
// Version  : v3.0                                                            //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// INCLUDE                                                                    //
// ========================================================================== //

// Standard Template Library
# include <iostream>

# include "surftri.hpp"

// ========================================================================== //
// MAIN                                                                       //
// ========================================================================== //
int main(
        void
        ) {

    // ========================================================================== //
    // int main(                                                                  //
    //     void)                                                                  //
    //                                                                            //
    // Demo with example of usage for Class_SurfTri variables. Main.              //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //


    // Counters
    // none

    // Output message ------------------------------------------------------- //
    cout << "===================== Class_SurfTri DEMO ===================== " << endl;
    // Scope variables ------------------------------------------------------ //
    // none

    // Run demo ------------------------------------------------------------- //
    Demo_CleaningTools(); 
    Demo_CleaningTools2();
    Demo_GenerationTools_2D();
    Demo_GenerationTools_3D();
    Demo_CleaningTools3();

    return 0;
}

// -------------------------------------------------------------------------- //
void Demo_CleaningTools(
        void
        ) {

    // ========================================================================== //
    // void Demo_CleaningTools(                                                   //
    //     void)                                                                  //
    //                                                                            //
    // Demo for Class_SurfTri cleaning tools.                                     //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    Class_SurfTri           Tri;
    ivector2D		E, S2E;

    // Counters
    // none

    // ========================================================================== //
    // OUTPUT MESSAGE                                                             //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << "============= Class_SurfTri: cleaning tools demo ==============" << endl;
    }

    // ========================================================================== //
    // LOAD TASSELATION FROM DGF FILE                                             //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Loading surface tasselation" << endl << endl;

        // Import from dgf file ------------------------------------------------- //
        Tri.Import_dgf("./datasets/naca0012.dgf");

        // Output stats --------------------------------------------------------- //
        // Tri.Stats(cout);

        // Export to vtk format ------------------------------------------------- //
        Tri.Export_vtu("tri_0.vtu");

    }

    // ========================================================================== //
    // REMOVE DOUBLE VERTICES                                                     //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Removing duplicated vertices" << endl << endl;


        // Output stats --------------------------------------------------------- //
        Tri.RemoveDoubleVertex();

        cout << " test" << endl ;
        // Output stats --------------------------------------------------------- //
        // Tri.Stats(cout);


        // Export to vtk format ------------------------------------------------- //
        Tri.Export_vtu("tri_1.vtu");
    }

    // ========================================================================== //
    // REMOVE ISOLATED SIMPLICIES                                                 //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Removing isolated simplicies" << endl << endl;

        // Output stats --------------------------------------------------------- //
        Tri.RemoveIsolatedSimplex();

        // Output stats --------------------------------------------------------- //
        // Tri.Stats(cout);

        // Export to vtk format ------------------------------------------------- //
        Tri.Export_vtu("tri_2.vtu");
    }

    // ========================================================================== //
    // REMOVE ISOLATED VERTICES                                                   //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Removing isolated vertices" << endl << endl;

        // Output stats --------------------------------------------------------- //
        Tri.RemoveIsolatedVertex();

        // Output stats --------------------------------------------------------- //
        // Tri.Stats(cout);

        // Export to vtk format ------------------------------------------------- //
        Tri.Export_vtu("tri_3.vtu");
    }

    // ========================================================================== //
    // BUILD ADJACENCY                                                            //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        ivector1D           list;

        // Output message ------------------------------------------------------- //
        cout << endl << " - Building adjacency matrix" << endl << endl;

        // Output stats --------------------------------------------------------- //
        Tri.BuildAdjacency();

        // Output stats --------------------------------------------------------- //
        // Tri.Stats(cout);
        for (int i = 0; i < Tri.nSimplex; i++) {
            cout << Tri.Adjacency[i] << endl;
        } //next i

        // Output message ------------------------------------------------------- //
        cout << endl << " - Updating adjacency matrix" << endl << endl;

        // Update adjacencies --------------------------------------------------- //    
        list.resize(2, 0);
        list[1] = 1;
        Tri.UpdateAdjacency(list);
        for (int i = 0; i < Tri.nSimplex; i++) {
            cout << Tri.Adjacency[i] << endl;
        } //next i
    }

    // ========================================================================== //
    // REMOVE DOUBLE SIMPLICIES                                                   //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        ivector1D           list;

        // Output message ------------------------------------------------------- //
        cout << endl << " - Removing double simplicies" << endl << endl;

        // Output stats --------------------------------------------------------- //
        list = Tri.FindDoubleSimplex();
        Tri.RemoveSimplex(list);

        // Output stats --------------------------------------------------------- //
        // Tri.Stats(cout);
        for (int i = 0; i < Tri.nSimplex; i++) {
            cout << Tri.Adjacency[i] << endl;
        } //next i

        // Export to vtk format ------------------------------------------------- //
        Tri.Export_vtu("tri_4.vtu");

    }


    // ========================================================================== //
    // CLOSING MESSAGE                                                            //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << "======================== DEMO: done!! =========================" << endl;
    }

    return; };

// -------------------------------------------------------------------------- //
void Demo_CleaningTools2(
        void
        ) {

    // ========================================================================== //
    // void Demo_CleaningTools2(                                                  //
    //     void)                                                                  //
    //                                                                            //
    // Demo for Class_SurfTri cleaning tools.                                     //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    Class_SurfTri           Tri;

    // Counters
    // none

    // ========================================================================== //
    // OUTPUT MESSAGE                                                             //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << "============ Class_SurfTri: cleaning tools demo 2 =============" << endl;
    }

    // ========================================================================== //
    // LOAD TASSELATION FROM DGF FILE                                             //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Loading surface tasselation" << endl << endl;

        // Import from dgf file ------------------------------------------------- //
        Tri.Import_stl("./datasets/buddha.stl",true);

        // Output stats --------------------------------------------------------- //
        //     Tri.Stats(cout);

        // Export to vtk format ------------------------------------------------- //
        //     Tri.Export_vtu("tri_0.vtu");

    }

    // ========================================================================== //
    // REMOVE DUPLICATED VERTICES                                                 //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        Class_SurfTri                                Tri2 = Tri;
        chrono::time_point<chrono::system_clock>     start, end;
        int                                          elapsed_seconds;

        // Output message ------------------------------------------------------- //
        cout << endl << " - Removing duplicated vertices" << endl << endl;

        // Output stats --------------------------------------------------------- //
        start = std::chrono::system_clock::now();
        Tri2.RemoveDoubleVertex();
        Tri2.ResizeVertex();
        end = chrono::system_clock::now();
        elapsed_seconds = chrono::duration_cast<chrono::milliseconds>(end-start).count();
        cout << "elapsed time: " << elapsed_seconds << " ms" << endl;
        cout << "# of vertices:" << endl;
        cout << " before cleaning: " << Tri.nVertex << ", after cleaning: " << Tri2.nVertex << endl;
        cout << "# of simplicies:" << endl;
        cout << " before cleaning: " << Tri.nSimplex << ", after cleaning: " << Tri2.nSimplex << endl;;
        cout << " double vertices found after cleaning: " << Tri2.CountDoubleVertex() << endl;

        // Measure error w.r.t original surface --------------------------------- //
        double err = 0.0;
        int    m;
        for (int T = 0; T < Tri.nSimplex; T++) {
            m = Tri.Simplex[T].size();
            for (int i = 0; i < m; i++) {
                err += norm_2(Tri.Vertex[Tri.Simplex[T][i]] - Tri2.Vertex[Tri2.Simplex[T][i]]);
            } //next i
        } //next T
        cout << "error: " << err << endl;

        // Output stats --------------------------------------------------------- //
        //     Tri.Stats(cout);

        // Export to vtk format ------------------------------------------------- //
        //     Tri.Export_vtu("tri_1.vtu");
        //      Tri2.Export_vtu("tri_1b.vtu");
    }

    // ========================================================================== //
    // BUILD ADJACENCIES                                                          //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        Class_SurfTri                                Tri2 = Tri;
        chrono::time_point<chrono::system_clock>     start, end;
        int                                          elapsed_seconds;

        // Output message ------------------------------------------------------- //
        cout << endl << " - Building adjacencies" << endl << endl;

        // Output stats --------------------------------------------------------- //
        start = std::chrono::system_clock::now();
        Tri2.BuildAdjacency();
        end = chrono::system_clock::now();
        elapsed_seconds = chrono::duration_cast<chrono::milliseconds>(end-start).count();
        cout << "elapsed time: " << elapsed_seconds << " ms" << endl;

    }

    // ========================================================================== //
    // CLOSING MESSAGE                                                            //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << "======================== DEMO: done!! =========================" << endl;
    }

    return; };

// -------------------------------------------------------------------------- //
void Demo_GenerationTools_3D(
        void
        ) {

    // ========================================================================== //
    // void Demo_GenerationTools_3D(                                              //
    //     void)                                                                  //
    //                                                                            //
    // Demo for generation tools.                                                 //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    bool                    stl_type = false;
    string                  stl_name = "./datasets/cube.stl";
    Class_SurfTri           Tri, Tri2;
    Class_SurfTri           Edges;


    // ========================================================================== //
    // OUTPUT MESSAGE                                                             //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << "============ Class_SurfTri: generation tools demo (3D) ========" << endl;
    }

    // ========================================================================== //
    // LOAD TASSELATION FROM DGF FILE                                             //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Loading 3D surface tasselation" << endl << endl;

        // Import from dgf file ------------------------------------------------- //
        Tri.Import_stl(stl_name, stl_type);

        // Export to vtk format ------------------------------------------------- //
        Tri.Export_vtu("geom.vtu");

    }

    // ========================================================================== //
    // GENERATE EDGES                                                             //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Cleaning input triangulation" << endl << endl;

        // Clean input triangulation -------------------------------------------- //
        Tri.RemoveDoubleVertex();
        Tri.ResizeVertex();
        Tri.BuildAdjacency();

        // Build edges ---------------------------------------------------------- //

        // Build edge-vertex and simplex-edge connectivities
        Tri.BuildEdges();

        // Copy edge network into a surf_tri class
        Edges.AddVertices(Tri.Vertex);
        Edges.AddSimplicies(Tri.Edge);

        // Export edge network to .vtu file
        Edges.Export_vtu("edges.vtu");

    }

    // ========================================================================== //
    // GENERATE EDES NORMALS                                                      //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Generating edges normals" << endl << endl;

        // Generate edges normals ----------------------------------------------- //
        Tri.GenerateENormals();

        // Export results ------------------------------------------------------- //
        Edges.ExportVCData_vtu("Enormals.vtu", "n", Tri.ENormal);

    }

    // ========================================================================== //
    // GENERATE VERTEX NORMALS                                                    //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Generating vertex normals" << endl << endl;

        // Generate edges normals ----------------------------------------------- //
        Tri.GenerateVNormals();

        // Export results ------------------------------------------------------- //
        Tri.ExportVPData_vtu("Vnormals.vtu", "n", Tri.VNormal);

    }

    // ========================================================================== //
    // TEST ASSIGNAMENT OPERATORS                                                 //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Testing assignament operators" << endl << endl;

        // Makes a copy --------------------------------------------------------- //
        Tri2 = Tri;

        // Export edge and vertex normals --------------------------------------- //
        Edges.ExportVCData_vtu("Enormals_copy.vtu", "n", Tri2.ENormal);
        Tri2.ExportVPData_vtu("Vnormals_copy.vtu", "n", Tri2.VNormal);

    }

    // ========================================================================== //
    // TEST EDGE CONSTRUCTION WITH EXTERNAL VERTEX LIST                           //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Testing edge/vertex normals construction with external vertex list" << endl << endl;

        // Clean data structure ------------------------------------------------- //
        Tri2.DestroyVertex();
        Tri2.DestroyEdge();
        Tri2.DestroySimplex2Edge();
        Tri2.DestroyVNormal();
        Tri2.DestroyENormal();

        // Re-compute vertex normals -------------------------------------------- //
        Tri2.GenerateVNormals(Tri.Vertex);
        Tri2.AddVertices(Tri.Vertex);

        // Export re-generated edge and vertex normals -------------------------- //
        Edges.ExportVCData_vtu("Enormals_copy2.vtu", "n", Tri2.ENormal);
        Tri2.ExportVPData_vtu("Vnormals_copy2.vtu", "n", Tri2.VNormal);

    }


    // ========================================================================== //
    // CLOSING MESSAGE                                                            //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << "======================== DEMO: done!! =========================" << endl;
    }

    return; };

// -------------------------------------------------------------------------- //
void Demo_GenerationTools_2D(
        void
        ) {

    // ========================================================================== //
    // void Demo_GenerationTools_2D(                                              //
    //     void)                                                                  //
    //                                                                            //
    // Demo for generation tools.                                                 //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    string                  dgf_name = "./datasets/naca0012.dgf";
    Class_SurfTri           Tri, Tri2;
    Class_SurfTri           Edges;


    // ========================================================================== //
    // OUTPUT MESSAGE                                                             //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << "============ Class_SurfTri: generation tools demo (2D) ========" << endl;
    }

    // ========================================================================== //
    // LOAD TASSELATION FROM DGF FILE                                             //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Loading 2D surface tasselation" << endl << endl;

        // Import from dgf file ------------------------------------------------- //
        Tri.Import_dgf(dgf_name);

        // Export to vtk format ------------------------------------------------- //
        Tri.Export_vtu("geom.vtu");

    }

    // ========================================================================== //
    // GENERATE EDGES                                                             //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Cleaning input triangulation" << endl << endl;

        // Clean input triangulation -------------------------------------------- //
        Tri.RemoveDoubleVertex();
        Tri.ResizeVertex();
        Tri.BuildAdjacency();
        Tri.GenerateNormals();
        for (int i = 0; i < Tri.nSimplex; ++i) {
            swap(Tri.Normal[i][0], Tri.Normal[i][1]);
            Tri.Normal[i][1] = -Tri.Normal[i][1];
        } //next i

        // Build edges ---------------------------------------------------------- //

        // Build edge-vertex and simplex-edge connectivities
        Tri.BuildEdges();

        // Copy edge network into a surf_tri class
        Edges.AddVertices(Tri.Vertex);
        Edges.AddSimplicies(Tri.Edge);

        // Export edge network to .vtu file
        Edges.Export_vtu("edges.vtu");

    }

    // ========================================================================== //
    // GENERATE EDES NORMALS                                                      //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Generating edges normals" << endl << endl;

        // Generate edges normals ----------------------------------------------- //
        Tri.GenerateENormals();

        // Export results ------------------------------------------------------- //
        Edges.ExportVCData_vtu("Enormals.vtu", "n", Tri.ENormal);

    }

    // ========================================================================== //
    // GENERATE VERTEX NORMALS                                                    //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Generating vertex normals" << endl << endl;

        // Generate edges normals ----------------------------------------------- //
        Tri.GenerateVNormals();

        // Export results ------------------------------------------------------- //
        Tri.ExportVPData_vtu("Vnormals.vtu", "n", Tri.VNormal);

    }

    // ========================================================================== //
    // TEST ASSIGNAMENT OPERATORS                                                 //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Testing assignament operators" << endl << endl;

        // Makes a copy --------------------------------------------------------- //
        Tri2 = Tri;

        // Export edge and vertex normals --------------------------------------- //
        Edges.ExportVCData_vtu("Enormals_copy.vtu", "n", Tri2.ENormal);
        Tri2.ExportVPData_vtu("Vnormals_copy.vtu", "n", Tri2.VNormal);

    }

    // ========================================================================== //
    // TEST EDGE CONSTRUCTION WITH EXTERNAL VERTEX LIST                           //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Testing edge/vertex normals construction with external vertex list" << endl << endl;

        // Clean data structure ------------------------------------------------- //
        Tri2.DestroyVertex();
        Tri2.DestroyEdge();
        Tri2.DestroySimplex2Edge();
        Tri2.DestroyVNormal();
        Tri2.DestroyENormal();

        // Re-compute vertex normals -------------------------------------------- //
        Tri2.GenerateVNormals(Tri.Vertex);
        Tri2.AddVertices(Tri.Vertex);

        // Export re-generated edge and vertex normals -------------------------- //
        Edges.ExportVCData_vtu("Enormals_copy2.vtu", "n", Tri2.ENormal);
        Tri2.ExportVPData_vtu("Vnormals_copy2.vtu", "n", Tri2.VNormal);

    }


    // ========================================================================== //
    // CLOSING MESSAGE                                                            //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << "======================== DEMO: done!! =========================" << endl;
    }

    return; };

// -------------------------------------------------------------------------- //
void Demo_CleaningTools3(
        void
        ) {

    // ========================================================================== //
    // void Demo_CleaningTools3(                                                  //
    //     void)                                                                  //
    //                                                                            //
    // Demo for Class_SurfTri cleaning tools. check on TrueDoubleSimplex and      //
    // Ring with safe failure are performed
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    Class_SurfTri           Tri;
    ivector2D		E, S2E;
    std::string selection, tag;
    // Counters
    // none

    // ========================================================================== //
    // OUTPUT MESSAGE                                                             //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << "============= Class_SurfTri: cleaning tools 3 demo ==============" << endl;
    }

    // ========================================================================== //
    // LOAD TASSELATION FROM STL FILE                                             //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none


        // Output message ------------------------------------------------------- //
        cout << " - Reading surface tasselation :" << endl << endl;
        cout << " Enter the name of your geometry .stl/STL (absolute path needed)  :";
        cin >> selection;


        tag = selection.substr(selection.size()-4,4);
        if(tag !=".stl" && tag != ".STL") {cout<<"Geometry not supported in this demo. Now exiting...."<<endl; exit(1);}

        cout<<"Going to read "<<selection<<" with tag"<< tag <<endl;
        // Import from dgf file ------------------------------------------------- //
        Tri.Import_stl(selection, tag==".STL");

        cout<<"number of vertices  "<<Tri.nVertex<<endl;  
        cout<<"number of simplicies  "<<Tri.nSimplex<<endl;  
        // Output stats --------------------------------------------------------- //

        Tri.SetTolerance();				

    }

    // ========================================================================== //
    // REMOVE DOUBLE VERTICES                                                     //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Finding duplicated vertices" << endl;

        // Output stats --------------------------------------------------------- //
        ivector1D list = Tri.FindDoubleVertex();  
        Tri.RemoveDoubleVertex();
        cout << endl << " - Going to remove  " << list.size() <<"   vertices"<<endl<<endl;				
        Tri.BuildAdjacency();	
    }

    // ========================================================================== //
    // REMOVE 0 Area SIMPLICIES                                                 //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Finding 0 area simplicies" << endl;

        ivector1D list;
        // Output stats --------------------------------------------------------- //
        list = Tri.Find0AreaSimplex();
        cout << endl << " - Going to remove  " << list.size() << " 0 area simplicies"<<endl<<endl;				
        Tri.RemoveSimplex(list);
        Tri.UpdateAdjacency(list);
    }

    // ========================================================================== //
    // REMOVE ISOLATED SIMPLICIES                                                 //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Finding isolated simplicies" << endl;

        ivector1D list;
        // Output stats --------------------------------------------------------- //
        list = Tri.FindIsolatedSimplex();
        cout << endl << " - Going to remove  " << list.size() << "  simplicies"<<endl<<endl;				
        Tri.RemoveSimplex(list);
        Tri.UpdateAdjacency(list);
    }
    // ========================================================================== //
    // REMOVE ISOLATED VERTICES                                                   //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Finding isolated vertices" << endl;

        ivector1D list;
        // Output stats --------------------------------------------------------- //
        list = Tri.FindIsolatedVertex();
        cout << endl << " - Going to remove  " << list.size() << "  vertices"<<endl<<endl;				
        Tri.RemoveVertex(list);
        Tri.UpdateAdjacency(list);
    }


    // ========================================================================== //
    // REMOVE DOUBLE SIMPLICIES                                                   //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        ivector1D           list;

        // Output message ------------------------------------------------------- //
        cout << endl << " - Removing true double simplicies" << endl << endl;

        // Output stats --------------------------------------------------------- //
        list = Tri.FindTrueDoubleSimplex();
        cout << endl << " - Going to remove  " << list.size() <<"  simplicies"<<endl<<endl;	
        Tri.RemoveSimplex(list);
        Tri.UpdateAdjacency(list);
    }

    Tri.FixNodeNumb();
    Tri.ResizeVertex();
    Tri.ResizeSimplex();
    Tri.ResizeNormal();
    Tri.ResizeAdjacency();

    Tri.Stats(cout);
    Tri.Export_vtu("cleanedgrid.vtu");

    // ========================================================================== //
    // CLOSING MESSAGE                                                            //
    // ========================================================================== //
    {
        //     // Scope variables ------------------------------------------------------ //
        //     // none
        // 
        //     ivector2D edges, edgesadj;
        //     Tri.BuildEdges(edges,edgesadj);
        // 
        //     int counter = 0;
        //     cout<<"Check every 1-ring in tassellation vertices...."<<endl;
        // 
        //     bvector1D checkV(Tri.nVertex, false);
        // 				
        //    for(int T=0; T<Tri.nSimplex; ++T)
        //    {
        //     for(int j=0; j<Tri.Simplex[0].size(); ++j)
        //     {
        //      bool check,isRing=true;
        //      if(!checkV[Tri.Simplex[T][j]])
        // 	{										
        //          ivector1D list = Tri.Ring_1(T,j,check, isRing);
        // 	 if(!isRing) {	cout<<"Failed Ring Computation for vertex "<<Tri.Simplex[T][j]<<" ring size was"<<list.size()<<endl;
        // 			counter++;								     
        // 	             }	
        //      	checkV[Tri.Simplex[T][j]] = true;
        //   	}	
        //     }
        //   }

        //  cout<<"Total failed Rings caught:   "<<counter; 	
        //  if(counter>0) {cout<<"   Check your triangulation."<<endl;}
        // Output message ------------------------------------------------------- //
        cout << "======================== DEMO: done!! =========================" << endl;
    }

    return; };

