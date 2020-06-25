#include "utilities/utility_functions.h"

using namespace utility_functions;

template<class T> void load_tree(T& tree, cli_parameters &cli);
template<class T> void morse_simplification(T& tree, cli_parameters &cli);

int main(int , char** )
{
    cli_parameters cli;
    cli.mesh_path = "../data/devil_0.tri";
    cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = 10;
	cli.app_debug = OUTPUT;
	cli.persistence = 0.8;
		
    cerr<<"[OBJECTIVE] this unit-test generates a quadtrees based on the PR-T tree criterion. "
	    <<"then, it saves the index and the mesh in VTK format for visualization purposes and finally "
		<<"it computes the Morse gradient vector and simplifies it, outputting the initial Morse IG and the simplified one in VTK format."<<endl;
		
    cerr<<"[NOTA] all the generated files are saved in the 'data' folder"<<endl;
	
	setParameters(cli);
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PR-T tree"<<endl;
    load_tree(ptree,cli);
    morse_simplification(ptree,cli);    

    return (EXIT_SUCCESS);
}

template<class T> void load_tree(T& tree, cli_parameters &cli)
{
    Timer time;
    
    
    if (!Reader::read_mesh(tree.get_mesh(), cli.mesh_path))
    {
        cout << "[ERROR] Loading mesh file. Execution Stopped." << endl;
        return;
    }

    cerr << "[INFO] Generating a Terrain tree: " << cli.v_per_leaf << " " << cli.t_per_leaf << " " << cli.crit_type << " ";
    time.start();
    tree.build_tree();
    time.stop();
    time.print_elapsed_time("[TIME] Building: ");

    cli.rever_to_original = Reader::read_noised_field_value(tree.get_mesh(),get_path_without_file_extension(cli.mesh_path),cli.original_vertex_fields);

    cli.original_vertex_indices.assign(tree.get_mesh().get_vertices_num(),-1);
    cli.original_triangle_indices.assign(tree.get_mesh().get_triangles_num(),-1);

    time.start();
    Reindexer reindexer = Reindexer();
    reindexer.reindex_tree_and_mesh(tree,true,cli.original_vertex_indices,
                                    true,cli.original_triangle_indices);
    time.stop();
    time.print_elapsed_time("Index and Mesh Reindexing ");
}

template<class T> void morse_simplification(T& tree, cli_parameters &cli)
{
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);

    //CALCOLO IL FORMAN GRADIENT VECTOR
    Forman_Gradient forman_gradient = Forman_Gradient(tree.get_mesh().get_triangles_num());
    Forman_Gradient_Computation gradient_computation = Forman_Gradient_Computation();

    Timer time = Timer();

    /// ---- FORMAN GRADIENT COMPUTATION --- ///
    cout<<"[NOTA] Compute the gradient field"<<endl;
    time.start();
    gradient_computation.compute_gradient_vector(forman_gradient,tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] computing gradient vector field ");

    /// ---- MORPHOLOGICAL SIMPLIFICATION --- ///        
    {
        Forman_Gradient_Simplifier forman_simplifier;

        ///
        /// firstly we extract the MIG
        ///
        cout<<"--- Morse Incidence Graph BEFORE simplification ---"<<endl;
        forman_simplifier.get_incidence_graph().init(); /// init again the base of the MIG
        forman_simplifier.extract_incidence_graph(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),OUTPUT,cli.cache_size); /// we force to keep the MIG structure


        Writer_Morse::write_incidence_graph_VTK(out.str(),"mig", cli.v_per_leaf, forman_simplifier.get_incidence_graph(),tree.get_mesh(),
                                          cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original); /// and we save it

        forman_simplifier.reset_stats();
        forman_simplifier.reset_output_structures(tree.get_mesh());
        ///
        /// then we execute effectively the topological simplification
        ///
        /// we have chosen a fully local simplification. i.e. the MIG is computed locally
        {
            cout<<"[LOCALLY] Simplify the forman gradient vector."<<endl;
            time.start();
            forman_simplifier.exec_local_topological_simplification(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),
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
        forman_simplifier.reset_stats();

        Writer_Morse::write_incidence_graph_VTK(out.str(),"simplified_mig", cli.v_per_leaf, forman_simplifier.get_incidence_graph(),tree.get_mesh(),
                                          cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);

        forman_simplifier.reset_output_structures(tree.get_mesh());
        forman_simplifier.reset_timer_variables();
    }
}
