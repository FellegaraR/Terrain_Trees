#include "utilities/utility_functions.h"

using namespace utility_functions;

void load_tree(PRT_Tree& tree, cli_parameters &cli);

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
    cerr<<"[GENERATION] PR-T tree"<<endl;
    if(strcmp(argv[3],"-q")==0)
      cli.QEM_based=true;
    else
      cli.QEM_based=false;
    load_tree(ptree,cli);

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
    
    tree.init_leaves_list(tree.get_root()); // NOTA: enables the parallel OpenMP-based processing of the leaf block

    Statistics stats;
    stats.get_index_statistics(tree,cli.reindex);

    Contraction_Simplifier simplifier;
    simplifier.simplify_parallel(tree,tree.get_mesh(),cli);

    cout<<"number of remaining triangles: "<<tree.get_mesh().get_triangles_num()<<endl;

    
    cout<<output_name<<endl;
    Writer::write_mesh_VTK(output_name,tree.get_mesh());  

}
