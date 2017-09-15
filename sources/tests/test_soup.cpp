#include "utilities/utility_functions.h"

using namespace utility_functions;

template<class T> void load_tree(T& tree, cli_parameters &cli);

int main(int , char** )
{
	cli_parameters cli;
	cli.mesh_path = "../data/simple_terrain.soup";
	
	cerr<<"[OBJECTIVE] this unit-test generates a quadtree from a soup of triangles instead of it classical indexed representation. "
	    <<" It saves, in VTK format outputs "
	    <<" the index and the mesh for visualization purposes and finally "
		<<" it computes the spatial index statistics."<<endl;
		
    cerr<<"[NOTA] all the generated files are saved in the 'data' folder"<<endl;
	
	cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = 2;
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PR-T tree"<<endl;
    load_tree(ptree,cli);

    return (EXIT_SUCCESS);
}

template<class T> void load_tree(T& tree, cli_parameters &cli)
{
    Timer time;
    Soup soup;

    //Legge l'input
    if(get_file_extension(cli.mesh_path) == "soup")
    {
        /// triangle mesh as soup
        if (!Reader::read_soup(soup, cli.mesh_path))
        {
            cerr << "[ERROR] Loading soup file. Execution Stopped." << endl;
            return;
        }
    }

    stringstream base_info;
    base_info << cli.v_per_leaf << " " << cli.t_per_leaf << " " << cli.crit_type << " ";
    stringstream base;
    base << get_path_without_file_extension(cli.mesh_path);

    stringstream tree_info;
    tree_info << base_info.str() << "[TIME] Building ";
    time.start();
    tree.build_tree(soup);
    soup.clear();
    time.stop();
    time.print_elapsed_time(tree_info.str());

    stringstream out2;
    out2 << base.str();
    out2 << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type << "_v_" << cli.v_per_leaf << "_tree.vtk";

    Writer::write_tree_VTK(out2.str(),tree.get_root(),tree.get_subdivision(),tree.get_mesh());
    Writer::write_mesh_VTK(base.str(),tree.get_mesh());        

    time.start();
    Reindexer reindexer = Reindexer();
    reindexer.reindex_tree_and_mesh(tree,false,cli.original_vertex_indices,
                                    false,cli.original_triangle_indices);
    time.stop();
    time.print_elapsed_time("[TIME] Index and Mesh Reindexing ");
    
    Statistics stats;
    stats.get_index_statistics(tree,cli.reindex);
}
