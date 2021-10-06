#include "utilities/utility_functions.h"

using namespace utility_functions;

template<class T> void load_tree(T& tree, cli_parameters &cli);
template<class T> void multi_morse_terrain_analysis(T& tree, cli_parameters &cli);
void multi_morse_terrain_analysis(PMRT_Tree &tree, cli_parameters &cli);

int main(int , char** )
{
    cli_parameters cli;
    cli.mesh_path = "../data/eggs64z4.tri";
    cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = 10;
	cli.app_debug = OUTPUT;	
	
    cerr<<"[OBJECTIVE] this unit-test generates a quadtrees based on the PR-T tree and PMR-T tree criteria. "
        <<"It computes first the MultiVariate Morse gradient vector and finally it extracts the critical clusters, outputting them in a VTK file."<<endl;
    cerr<<"[NOTA] all the generated files are saved in the 'data' folder"<<endl;
	
	setParameters(cli);
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PR-T tree"<<endl;
    load_tree(ptree,cli);
    multi_morse_terrain_analysis(ptree,cli);

    cli.mesh_path = "../data/eggs64z4_extra_fields.off";
    cli.v_per_leaf = -1;
    cli.t_per_leaf = 20;
    cli.crit_type = "pmr";
    cerr<<"[GENERATION] PMR-T tree"<<endl;
    PMRT_Tree rttree = PMRT_Tree(cli.t_per_leaf,cli.division_type);
    load_tree(rttree,cli);
    multi_morse_terrain_analysis(rttree,cli);

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

    time.start();
    Reindexer reindexer = Reindexer();
    reindexer.reindex_tree_and_mesh(tree,false,cli.original_vertex_indices,
                                    false,cli.original_triangle_indices);
    time.stop();
    time.print_elapsed_time("[TIME] Index and Mesh Reindexing ");
}

template<class T> void multi_morse_terrain_analysis (T& tree, cli_parameters &cli)
{
    Timer time = Timer();

    C_Curvature ccurvature = C_Curvature(MEAN,true);
    time.start();
    ccurvature.compute(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] mean-cCurvature computation: ");

    Writer::write_mesh(get_path_without_file_extension(cli.mesh_path),"extra_fields",tree.get_mesh(),true);

    Forman_Gradient forman_gradient = Forman_Gradient(tree.get_mesh().get_triangles_num());
    Forman_Gradient_Computation gradient_computation = Forman_Gradient_Computation();

    cout<<"[NOTA] Compute the gradient field"<<endl;
    time.start();
    gradient_computation.compute_gradient_vector(forman_gradient,tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] computing gradient vector field ");

    /// this extraction must be executed at the end..
    /// as it invalidates the critical d.s.
    Forman_Gradient_Features_Extractor features_extractor;
    features_extractor.extract_critical_clusters(tree.get_root(),tree.get_mesh(),forman_gradient,gradient_computation.get_critical_simplices(),
                                                 tree.get_subdivision(), OUTPUT,cli.cache_size,get_path_without_file_extension(cli.mesh_path));
}

void multi_morse_terrain_analysis(PMRT_Tree &tree, cli_parameters &cli)
{
    Timer time = Timer();

    C_Curvature ccurvature = C_Curvature(MEAN,true);
    time.start();
    ccurvature.compute(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] mean-cCurvature computation: ");

    Forman_Gradient forman_gradient = Forman_Gradient(tree.get_mesh().get_triangles_num());
    Forman_Gradient_Computation gradient_computation = Forman_Gradient_Computation();

    cout<<"[NOTA] Compute the gradient field"<<endl;
    time.start();
    gradient_computation.compute_gradient_vector(forman_gradient,tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] computing gradient vector field ");

    /// this extraction must be executed at the end..
    /// as it invalidates the critical d.s.
    Forman_Gradient_Features_Extractor features_extractor;
    features_extractor.extract_critical_clusters(tree.get_root(),tree.get_mesh(),forman_gradient,gradient_computation.get_critical_simplices(),
                                                 tree.get_subdivision(), OUTPUT,cli.cache_size,get_path_without_file_extension(cli.mesh_path));
}
