#include "utilities/utility_functions.h"

using namespace utility_functions;

template<class T> void load_tree(T& tree, cli_parameters &cli);

int main(int , char** )
{
	cli_parameters cli;
	cli.mesh_path = "../data/eggs64z4.tri";
	
	cerr<<"[OBJECTIVE] this unit-test generates three quadtrees based on the three criteria available. "
	    <<"then, it saves the index and the mesh in VTK format for visualization purposes and finally "
		<<"it computes the spatial index statistics."<<endl;
		
    cerr<<"[NOTA] all the generated files are saved in the 'data' folder"<<endl;
	
	cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = 10;
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PR-T tree"<<endl;
    load_tree(ptree,cli);
    
    cli.crit_type = "pm";
    cli.t_per_leaf = 20;
    PMT_Tree pttree = PMT_Tree(cli.v_per_leaf, cli.t_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PT_Tree"<<endl;
    load_tree(pttree,cli);
    
    cli.v_per_leaf = -1;
    cli.crit_type = "pmr";
    cerr<<"[GENERATION] PMR-T tree"<<endl;
    PMRT_Tree rttree = PMRT_Tree(cli.t_per_leaf,cli.division_type);
    load_tree(rttree,cli);    

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

    stringstream out;
    out << base.str() << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type;
    if (cli.crit_type == "pr")
       out << "_v_" << cli.v_per_leaf << "_.tree";
    else if (cli.crit_type == "pm")
       out << "_v_" << cli.v_per_leaf << "_t_" << cli.t_per_leaf << "_.tree";
    else if (cli.crit_type == "pmr")
       out << "_t_" << cli.t_per_leaf << "_.tree";
    Writer::write_tree(out.str(), tree.get_root(), tree.get_subdivision());

    stringstream out2;
    out2 << base.str();
    if (cli.crit_type == "pr")
       out2 << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type << "_v_" << cli.v_per_leaf << "_tree.vtk";
    else if (cli.crit_type == "pm")
       out2 << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type << "_v_" << cli.v_per_leaf << "_t_" << cli.t_per_leaf << "_tree.vtk";
    else if (cli.crit_type == "pmr")
       out2 << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type << "_t_" << cli.t_per_leaf << "_tree.vtk";

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
