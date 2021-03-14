#include "utilities/utility_functions.h"

using namespace utility_functions;

void load_tree_lite(PRT_Tree& tree, cli_parameters &cli);
void gradient_aware_simplification(PRT_Tree& tree, cli_parameters &cli);
void load_terrain(PRT_Tree& tree, cli_parameters &cli);

int main(int argc, char** argv)
{
	cli_parameters cli;
	cli.mesh_path = argv[1];
	cli.debug_mode = false;
	cerr<<"[OBJECTIVE] this unit-test generates a PR-quadtree on the input TIN dataset "
	    <<"then, it simplifies the triangle mesh with an edge contraction operator following a length criteria."<<endl;
		
    cerr<<"[NOTA] all the generated files are saved in the 'data' folder"<<endl;
	
	cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = atoi(argv[2]);

    if(atof(argv[4])==-1){
        cli.contract_all_edges = true;
    }
    
    cli.maximum_limit=atof(argv[4]);

    cli.num_of_threads =atoi(argv[5]);
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
        out << "_v_" << cli.v_per_leaf << "_.tree";
 
        cli.tree_path=out.str();

        if (!Reader::read_tree(tree, tree.get_root(), cli.tree_path))
        {
            cerr << "[ERROR] Loading .tree file." << endl;
            cerr << "[GENERATION] tree from triangle mesh" << endl;
            time.start();
            tree.build_tree();
            time.stop();
            time.print_elapsed_time(tree_info.str());
            Writer::write_tree(out.str(), tree.get_root(), tree.get_subdivision());
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


        time.start();
        Reindexer reindexer = Reindexer();
        reindexer.reindex_tree_and_mesh(tree,true,cli.original_vertex_indices,
                                        false,cli.original_triangle_indices);
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
    cerr << "[MEMORY] peak for computing gradient vector field: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

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

    tree.init_leaves_list(tree.get_root()); 

    cout<<"[NOTA]Border checking"<<endl;
    Border_Checker border_checker=Border_Checker();
    border_checker.compute_borders(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] Border Checking ");
    time.start();
    Gradient_Aware_Simplifier simplifier;
    simplifier.preprocess(tree,tree.get_mesh(),cli);
    time.stop();
    time.print_elapsed_time("[TIME] Preporcessing");


      time.start();
    simplifier.gradient_aware_simplify_parallel(tree,tree.get_mesh(),cli,forman_gradient);
    time.stop();
    time.print_elapsed_time("[TIME] Gradient-aware simplification ");

    cout<<"number of remaining triangles: "<<tree.get_mesh().get_triangles_num()<<endl;

    
    // cout<<output_name<<endl;
    // Writer::write_mesh_VTK(output_name,tree.get_mesh()); 
    // Writer::write_mesh(output_name,"grad",tree.get_mesh(),false); 

}

