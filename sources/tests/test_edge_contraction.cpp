#include "utilities/utility_functions.h"

using namespace utility_functions;

template<class T> void load_tree(T& tree, cli_parameters &cli);

int main(int , char** )
{
	cli_parameters cli;
	cli.mesh_path = "../data/cos_sum.tri";
	
	cerr<<"[OBJECTIVE] this unit-test generates a PR-quadtree on the input TIN dataset "
	    <<"then, it simplifies the triangle mesh with an edge contraction operator following a length criteria."<<endl;
		
    cerr<<"[NOTA] all the generated files are saved in the 'data' folder"<<endl;
	
	cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = 50;
    cli.maximum_length=0.45;
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PR-T tree"<<endl;
    load_tree(ptree,cli);

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

   cout<<"number of triangles: "<<tree.get_mesh().get_triangles_num()<<endl;
cout<<"number of vertices: "<<tree.get_mesh().get_vertices_num()<<endl;
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
   
   string output_name=base.str()+"_simplified_v_"+to_string(cli.v_per_leaf);
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
