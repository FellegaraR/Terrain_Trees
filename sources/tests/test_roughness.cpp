#include "utilities/utility_functions.h"

using namespace utility_functions;

template<class T> void load_tree(T& tree, cli_parameters &cli);
template<class T> void compute_roughness(T& tree, cli_parameters &cli);
void compute_roughness(PMRT_Tree& tree, cli_parameters &cli);

int main(int argc , char** argv)
{
	cli_parameters cli;
    cli.mesh_path = argv[1];
    
    cerr<<"[OBJECTIVE] this unit-test generates a quadtrees based on the PR-T tree criterion. "
	    <<"then, it computes the roughness of all the vertices in the mesh. "<<endl;
		
    cerr<<"[NOTA] all the generated files are saved in the 'data' folder"<<endl;
    
    cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = 20;    
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PR-T tree"<<endl;
    load_tree(ptree,cli);
    compute_roughness(ptree,cli);

    cli.v_per_leaf = -1;
    cli.t_per_leaf = 50;
    PMRT_Tree rtree = PMRT_Tree(cli.t_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PMR-T tree"<<endl;
    load_tree(rtree,cli);
    compute_roughness(rtree,cli);

    cli.v_per_leaf = 20;
    cli.t_per_leaf = 50;
    PMT_Tree pttree = PMT_Tree(cli.v_per_leaf,cli.t_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PT_Tree"<<endl;
    load_tree(pttree,cli);
    compute_roughness(pttree,cli);

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
    tree_info << base_info.str() << "[TIME] Building ";
    time.start();
    tree.build_tree();
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
         cerr << "[MEMORY] peak for Index and Mesh Reindexing: " <<
        to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
}  

template<class T> void compute_roughness(T& tree, cli_parameters &cli)
{
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);
    
    cout<<"[NOTA]Border checking"<<endl;
    Border_Checker border_checker=Border_Checker();
    border_checker.compute_borders(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    cout<<"[NOTA] compute roughness"<<endl;
    Timer time;
    Roughness roughness = Roughness(tree.get_mesh());
    time.start();
    roughness.compute(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] roughness computation: ");
     cerr << "[MEMORY] peak for computing Roughness: " <<
        to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
           roughness.store_result(tree.get_mesh());
    roughness.print_roughness_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);

//    Writer::write_mesh_roughness_VTK(out.str(),tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
   }

void compute_roughness(PMRT_Tree& tree, cli_parameters &cli)
{
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);

    Timer time;
     cout<<"[NOTA]Border checking"<<endl;
    Border_Checker border_checker=Border_Checker();
    border_checker.compute_borders(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    cout<<"[NOTA] compute roughness"<<endl;
    Roughness roughness = Roughness(tree.get_mesh());
    time.start();
    roughness.compute(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] roughness computation: ");
         cerr << "[MEMORY] peak for computing Roughness: " <<
        to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
               roughness.store_result(tree.get_mesh());
    roughness.print_roughness_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
  //  Writer::write_mesh_roughness_VTK(out.str(),tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
    }
