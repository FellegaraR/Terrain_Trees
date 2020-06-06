/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "utilities/utility_functions.h"

using namespace utility_functions;

template<class T> void load_tree(T& tree, cli_parameters &cli);
template<class T> void exec_queries(T& tree);

int main(int argc, char** argv )
{
        cerr<<"[OBJECTIVE] this unit-test tests the index based on different kv and kt.  "
	    <<"then, it extracts VT relation of the stored tree. ";
     cli_parameters cli;
    cli.mesh_path = argv[1];
    cli.division_type = QUAD;
     cli.v_per_leaf = atoi(argv[2]);    
    cli.t_per_leaf= atoi(argv[3]);
        string tree_type=argv[4];
      if(tree_type=="PR")
    {
    cli.crit_type = "pr";
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PR-T tree"<<endl;
    load_tree(ptree,cli);
    exec_queries(ptree);
    
 
    }
    else if(tree_type=="PMR"){
            cli.crit_type = "pmr";
    cli.v_per_leaf = -1;

    PMRT_Tree rtree = PMRT_Tree(cli.t_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PMR-T tree"<<endl;
    load_tree(rtree,cli);
   exec_queries(rtree);  

    }
    else if(tree_type=="PM"){
            cli.crit_type = "pm";
    PMT_Tree pttree = PMT_Tree(cli.v_per_leaf,cli.t_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PM-T tree"<<endl;
    load_tree(pttree,cli);
    exec_queries(pttree);  

    }
    
  
     

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

    stringstream tree_info;
    cout<< base_info.str()<<"vertices number: "<<tree.get_mesh().get_vertices_num()<<" triangles number: " <<tree.get_mesh().get_triangles_num()<< "[TIME] Building "<<endl;

    
     stringstream out; 
     out << base.str() << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type;
    if (cli.crit_type == "pr")
        out << "_v_" << cli.v_per_leaf << "_.tree";
    else if (cli.crit_type == "pm")
        out << "_v_" << cli.v_per_leaf << "_t_" << cli.t_per_leaf << "_.tree";
    else if (cli.crit_type == "pmr")
        out << "_t_" << cli.t_per_leaf << "_.tree";
    cli.tree_path=out.str();
    
    time.start();
    if (!Reader::read_tree(tree, tree.get_root(), cli.tree_path))
    {
        cerr << "[ERROR] Loading .tree file. Build a new tree." << endl;
            tree.build_tree();
    }
    time.stop();
      time.print_elapsed_time(tree_info.str());
   
    //Writer::write_tree(out.str(), tree.get_root(), tree.get_subdivision());
    
    
    // stringstream out2;
    //out2 << base.str();
    //out2 << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type << "_v_" << cli.v_per_leaf << "_tree.vtk";
    
    //Writer::write_tree_VTK(out2.str(),tree.get_root(),tree.get_subdivision(),tree.get_mesh());
    //Writer::write_mesh_VTK(base.str(),tree.get_mesh());        

    time.start();
    Reindexer reindexer = Reindexer();
    reindexer.reindex_tree_and_mesh(tree,false,cli.original_vertex_indices,
                                    false,cli.original_triangle_indices);
    time.stop();
    time.print_elapsed_time("[TIME] Index and Mesh Reindexing ");

       Statistics stats;
    stats.get_index_statistics(tree,cli.reindex);
}  
template<class T> void exec_queries(T& tree)
{
    Topological_Queries tq;
    tq.batched_VT(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_subdivision(),true);
    tq.batched_VT(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_subdivision(),true);
    tq.batched_VT(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_subdivision(),true);
}
