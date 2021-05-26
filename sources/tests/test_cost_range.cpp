#include "utilities/utility_functions.h"

using namespace utility_functions;

void load_tree_lite(PRT_Tree& tree, cli_parameters &cli);
void cost_hist(PRT_Tree& tree, cli_parameters &cli, int num_bins);
void load_terrain(PRT_Tree& tree, cli_parameters &cli);
void count_critical_simplices(PRT_Tree& tree, cli_parameters &cli, Forman_Gradient& forman_gradient);


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

    int num_bins;
    num_bins = atoi(argv[3]);
    

    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cli.app_debug=OUTPUT;
    cerr<<"[GENERATION] PR-T tree"<<endl;

   cost_hist(ptree,cli, num_bins);
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

void cost_hist(PRT_Tree& tree, cli_parameters &cli, itype num_bins){

    stringstream out;
    stringstream base;
    base << get_path_without_file_extension(cli.mesh_path);
    load_terrain(tree,cli);

    Timer time = Timer();

    load_tree_lite(tree,cli);

    time.start();
    tree.init_leaves_list(tree.get_root()); 
    time.stop();
    time.print_elapsed_time("[TIME] Initialize leave list: ");
   // cout<<"[NOTA]Border checking"<<endl;

    time.start();
    Gradient_Aware_Simplifier simplifier;
    //coord_type min,max;
    simplifier.error_range(tree,tree.get_mesh(),cli,num_bins);
    time.stop();
    time.print_elapsed_time("[TIME] Counting error distribution");





}