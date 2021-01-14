#include "utilities/utility_functions.h"

using namespace utility_functions;

void load_tree(PRT_Tree& tree, cli_parameters &cli);
void load_tree_lite(PRT_Tree& tree, cli_parameters &cli);
void gradient_aware_simplification(PRT_Tree& tree, cli_parameters &cli);
void load_terrain(PRT_Tree& tree, cli_parameters &cli);


int main(int argc, char** argv)
{
	cli_parameters cli;
	cli.mesh_path = argv[1];
	
	cerr<<"[OBJECTIVE] this unit-test generates a PR-quadtree on the input TIN dataset "
	    <<"then, it simplifies the triangle mesh with an edge contraction operator following a length criteria."<<endl;
		
    cerr<<"[NOTA] all the generated files are saved in the 'data' folder"<<endl;
	
	cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = atoi(argv[2]);
    cli.maximum_limit=atof(argv[4]);
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cli.app_debug=OUTPUT;
    cerr<<"[GENERATION] PR-T tree"<<endl;
    if(strcmp(argv[3],"-q")==0)
      cli.QEM_based=true;
    else
      cli.QEM_based=false;
    //load_tree(ptree,cli);
   gradient_aware_simplification(ptree,cli);
    return (EXIT_SUCCESS);
}

void load_tree(PRT_Tree& tree, cli_parameters &cli)
{
    Timer time;
    if (!Reader::read_mesh(tree.get_mesh(), cli.mesh_path))
    {
        cout << "[ERROR] Loading mesh file. Execution Stopped." << endl;
        return;
    }

    stringstream base_info;
    base_info << cli.v_per_leaf << " " << cli.t_per_leaf << " " << cli.crit_type << " ";
    stringstream base;
    base << get_path_without_file_extension(cli.mesh_path);

    cerr << "[INFO] Generating a Terrain tree: " << cli.v_per_leaf << " " << cli.t_per_leaf << " " << cli.crit_type << " ";
    stringstream tree_info;
    tree_info << base_info.str() << "[TIME] Building ";
    time.start();
    tree.build_tree();
    time.stop();
    time.print_elapsed_time(tree_info.str());

   cout<<"number of triangles: "<<tree.get_mesh().get_triangles_num()<<endl;
cout<<"number of vertices: "<<tree.get_mesh().get_vertices_num()<<endl;
    stringstream out;
    out << base.str() << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type;
       out << "_v_" << cli.v_per_leaf << "_.tree";

    Writer::write_tree(out.str(), tree.get_root(), tree.get_subdivision());

    stringstream out2;
    out2 << base.str();
       out2 << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type << "_v_" << cli.v_per_leaf << "_tree.vtk";

   string output_name=base.str()+"_simplified_v_"+to_string(cli.v_per_leaf);
   if(cli.QEM_based==true){

      output_name=output_name+"_qem";

   }
   else{
       output_name=output_name+"_length";
   }
     Writer::write_tree_VTK(out2.str(),tree.get_root(),tree.get_subdivision(),tree.get_mesh());
    // Writer::write_mesh_VTK(base.str(),tree.get_mesh());        

    time.start();
    Reindexer reindexer = Reindexer();
    reindexer.reindex_tree_and_mesh(tree,false,cli.original_vertex_indices,
                                    false,cli.original_triangle_indices);
    time.stop();
    time.print_elapsed_time("[TIME] Index and Mesh Reindexing ");
    
    Statistics stats;
    stats.get_index_statistics(tree,cli.reindex);

    Contraction_Simplifier simplifier;
    simplifier.simplify(tree,tree.get_mesh(),cli);

    cout<<"number of remaining triangles: "<<tree.get_mesh().get_triangles_num()<<endl;

    
    cout<<output_name<<endl;
    Writer::write_mesh_VTK(output_name,tree.get_mesh());  

}

void load_terrain(PRT_Tree& tree, cli_parameters &cli)
{
    if (!Reader::read_mesh(tree.get_mesh(), cli.mesh_path))
    {
        cout << "[ERROR] Loading mesh file. Execution Stopped." << endl;
        return;
    }

    cerr << "[MEMORY] peak for Indexing the terrain: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
}

 void load_tree_lite(PRT_Tree& tree, cli_parameters &cli)
{
    Timer time;

    stringstream base_info;
    base_info << cli.v_per_leaf << " " << cli.t_per_leaf << " " << cli.crit_type << " ";
    stringstream base;
    base << get_path_without_file_extension(cli.mesh_path);

    stringstream tree_info;
    tree_info << base_info.str();

    if (cli.isTreeFile)
    {
        cout << "tree path: " << cli.tree_path << endl;
        time.start();
        if (!Reader::read_tree(tree, tree.get_root(), cli.tree_path))
        {
            cerr << "[ERROR] Loading .tree file. Execution Stopped." << endl;
            return;
        }
        time.stop();
        tree_info << "[TIME] Loading tree from file ";
        time.print_elapsed_time(tree_info.str());
    }
    else
    {
        tree_info << "[TIME] Building ";

        stringstream out;
        out << base.str() << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type;
        if (cli.crit_type == "pr")
            out << "_v_" << cli.v_per_leaf << "_.tree";
        else if (cli.crit_type == "pm")
            out << "_v_" << cli.v_per_leaf << "_t_" << cli.t_per_leaf << "_.tree";
        else if (cli.crit_type == "pmr")
            out << "_t_" << cli.t_per_leaf << "_.tree";
        cli.tree_path=out.str();

        if (!Reader::read_tree(tree, tree.get_root(), cli.tree_path))
        {
            cerr << "[ERROR] Loading .tree file." << endl;
            cerr << "[GENERATION] tree from triangle mesh" << endl;
            time.start();
            tree.build_tree();
            time.stop();
            time.print_elapsed_time(tree_info.str());
         //   Writer::write_tree(out.str(), tree.get_root(), tree.get_subdivision());
        }
        else
            cout << "[NOTICE] Found corresponding .tree file. Loaded tree from file successfully"<<endl;



    }

    cerr << "[MEMORY] peak for encoding the Terrain tree: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;



    if(cli.reindex)
    {
        cerr<<"[REINDEXING] tree and triangle mesh"<<endl;


        //        if((cli.query_type == MORSE_ANALYSIS || cli.query_type == LOCAL_MORSE_SIMPLIFICATION || cli.query_type == GLOBAL_MORSE_SIMPLIFICATION)
        //                && cli.app_debug == OUTPUT)
        cli.original_vertex_indices.assign(tree.get_mesh().get_vertices_num(),-1);
        if(cli.query_type == MORSE_ANALYSIS && cli.app_debug == OUTPUT)
            cli.original_triangle_indices.assign(tree.get_mesh().get_triangles_num(),-1);

        time.start();
        Reindexer reindexer = Reindexer();
        reindexer.reindex_tree_and_mesh(tree,(cli.original_vertex_indices.size() == tree.get_mesh().get_vertices_num()),cli.original_vertex_indices,
                                        (cli.original_triangle_indices.size() == tree.get_mesh().get_triangles_num()),cli.original_triangle_indices);
        time.stop();
        time.print_elapsed_time("[TIME] Index and Mesh Reindexing ");
    }


    cerr << "[MEMORY] peak for reindexing the Terrain tree: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    if (cli.is_index)
    {
        Statistics stats;
        stats.get_index_statistics(tree,cli.reindex);
    }
}

void gradient_aware_simplification(PRT_Tree& tree, cli_parameters &cli){

    stringstream out;
    stringstream base;
    base << get_path_without_file_extension(cli.mesh_path);
    

    load_terrain(tree,cli);


    //CALCOLO IL FORMAN GRADIENT VECTOR
    Forman_Gradient forman_gradient = Forman_Gradient(tree.get_mesh().get_triangles_num());
    Forman_Gradient_Computation gradient_computation = Forman_Gradient_Computation();

    Timer time = Timer();

    time.start();
    gradient_computation.initial_filtering_IA(tree.get_mesh());
    time.stop();
    time.print_elapsed_time("[TIME] Initial filtering ");

    load_tree_lite(tree,cli);

    /// ---- FORMAN GRADIENT COMPUTATION --- ///
    gradient_computation.reset_filtering(tree.get_mesh(),cli.original_vertex_indices);

    cout<<"[NOTA] Computing the gradient field"<<endl;
    time.start();
    gradient_computation.compute_gradient_vector(forman_gradient,tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] computing gradient vector field ");

    stringstream out2;
    out2 << base.str();
    out2 << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type << "_v_" << cli.v_per_leaf << "_tree.vtk";

   string output_name=base.str()+"_simplified_v_"+to_string(cli.v_per_leaf);

   if(cli.QEM_based==true){

      output_name=output_name+"_qem";

   }
   else{
       output_name=output_name+"_length";
   }

    cout<<"[NOTA]Border checking"<<endl;
    time.start();
    Border_Checker border_checker=Border_Checker();
    border_checker.compute_borders(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
        time.stop();
    time.print_elapsed_time("[TIME] Border Checking ");

    Gradient_Aware_Simplifier simplifier;
      time.start();
    simplifier.gradient_aware_simplify(tree,tree.get_mesh(),cli,forman_gradient);
    time.stop();
    time.print_elapsed_time("[TIME] Gradient-aware simplification ");

    cout<<"number of remaining triangles: "<<tree.get_mesh().get_triangles_num()<<endl;

    
    cout<<output_name<<endl;
    Writer::write_mesh_VTK(output_name,tree.get_mesh()); 
    Writer::write_mesh(output_name,"simplified",tree.get_mesh(),false); 
    // Forman_Gradient_Features_Extractor features_extractor;   
    // features_extractor.extract_incidence_graph(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),cli.app_debug,cli.cache_size);
    // features_extractor.print_stats();
/*
    /// ---- MORPHOLOGICAL SIMPLIFICATION --- ///        
    {
        Forman_Gradient_Simplifier forman_simplifier;
        forman_simplifier.set_filtration_vec(gradient_computation.get_filtration());

        ///
        /// firstly we extract the MIG
        ///
        cout<<"--- Morse Incidence Graph BEFORE simplification ---"<<endl;
        forman_simplifier.get_incidence_graph().init(); /// init again the base of the MIG
        forman_simplifier.extract_incidence_graph(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),OUTPUT,cli.cache_size); /// we force to keep the MIG structure
        

        Writer_Morse::write_incidence_graph_VTK(out.str(),"mig", cli.v_per_leaf, forman_simplifier.get_incidence_graph(),tree.get_mesh(),
                                          cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original); /// and we save it

       // features_extractor.get_incidence_graph().print_stats(true);
        forman_simplifier.get_incidence_graph().print_stats(true);
        forman_simplifier.reset_stats();
        forman_simplifier.reset_output_structures(tree.get_mesh());
        ///
        /// then we execute effectively the topological simplification
        ///
        /// we have chosen a fully local simplification. i.e. the MIG is computed locally
        if (cli.query_type==LOCAL_MORSE_SIMPLIFICATION)
        {
            cout<<"[LOCALLY] Simplify the forman gradient vector."<<endl;
            time.start();
            forman_simplifier.exec_local_topological_simplification(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),
                                                                    cli.app_debug,cli.cache_size,cli.persistence);
            time.stop();
            time.print_elapsed_time("[TIME] simplify the gradient ");
        }
        else if(cli.query_type==GLOBAL_MORSE_SIMPLIFICATION){

             cout<<"[GLOBALLY] Simplify the forman gradient vector."<<endl;
                      time.start();
            /// otherwise we simplify the gradient computing first a global MIG and then simplifying it and the gradient
            /// default behaviour with alltime!
            forman_simplifier.exec_global_topological_simplification(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),
                                                                     cli.app_debug,cli.cache_size,cli.persistence);
            time.stop();
            time.print_elapsed_time("[TIME] simplify the gradient ");
        }
        forman_simplifier.print_simplification_stats();
    

        ///
        /// then we compute again and output the simplified mig
        ///
        cout<<"--- Morse Incidence Graph AFTER simplification ---"<<endl;
        forman_simplifier.reset_output_structures(tree.get_mesh());
        forman_simplifier.get_incidence_graph().init(); /// init again the MIG structures
        forman_simplifier.extract_incidence_graph(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),OUTPUT,cli.cache_size);
        forman_simplifier.get_incidence_graph().print_stats(true);
        forman_simplifier.reset_stats();

        if(cli.query_type==LOCAL_MORSE_SIMPLIFICATION)
        Writer_Morse::write_incidence_graph_VTK(out.str(),"simplified_mig_local", cli.v_per_leaf, forman_simplifier.get_incidence_graph(),tree.get_mesh(),
                                          cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);
        else if (cli.query_type==GLOBAL_MORSE_SIMPLIFICATION)
         Writer_Morse::write_incidence_graph_VTK(out.str(),"simplified_mig_global", cli.v_per_leaf, forman_simplifier.get_incidence_graph(),tree.get_mesh(),
                                          cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);

        forman_simplifier.reset_output_structures(tree.get_mesh());
        forman_simplifier.reset_timer_variables();
    }*/
}