#include "utilities/utility_functions.h"

using namespace utility_functions;

template<class T> void load_tree(T& tree, cli_parameters &cli);
template<class T> void compute_curvature(T& tree, cli_parameters &cli);
void compute_curvature(PMRT_Tree& tree, cli_parameters &cli);

int main(int , char** )
{
	cli_parameters cli;
    cli.mesh_path = "../data/devil_0.tri";
	
    cerr<<"[OBJECTIVE] this unit-test generates a quadtrees based on the PR-T tree criterion. "
	    <<"then, it saves the index and the mesh in VTK format for visualization purposes and finally "
		<<"it computes the concentrated, the meanC and GaussC curvatures, outputting in three VTK files the mesh with the curvature fields."<<endl;
		
    cerr<<"[NOTA] all the generated files are saved in the 'data' folder"<<endl;
	
	cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = 20;
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PR-T tree"<<endl;
    load_tree(ptree,cli);
    compute_curvature(ptree,cli);

    cli.v_per_leaf = -1;
    cli.t_per_leaf = 50;
    PMRT_Tree rtree = PMRT_Tree(cli.t_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PMR-T tree"<<endl;
    load_tree(rtree,cli);
    compute_curvature(rtree,cli);

    cli.v_per_leaf = 20;
    cli.t_per_leaf = 50;
    PMT_Tree pttree = PMT_Tree(cli.v_per_leaf,cli.t_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PT_Tree"<<endl;
    load_tree(pttree,cli);
    compute_curvature(pttree,cli);

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
}

template<class T> void compute_curvature(T& tree, cli_parameters &cli)
{
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);

    cout<<"[NOTA] compute concentrated-curvature"<<endl;
    Timer time;
    Concentrated_Curvature curvature = Concentrated_Curvature(CONCENTRATED);
    time.start();
    curvature.compute_borders(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] border computation: ");
    time.start();
    curvature.compute(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] concentrated-curvature computation: ");
    curvature.print_curvature_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
    Writer::write_mesh_curvature_VTK(out.str(),tree.get_mesh(),"concentrated",tree.get_mesh().get_vertex(1).get_fields_num()-1);

    cout<<"[NOTA] compute mean-cCurvature"<<endl;
    C_Curvature ccurvature = C_Curvature(MEAN,true);
    time.start();
    ccurvature.compute(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] mean-cCurvature computation: ");
    ccurvature.print_curvature_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
    Writer::write_mesh_curvature_VTK(out.str(),tree.get_mesh(),"mean_cCurvature_sumAll",tree.get_mesh().get_vertex(1).get_fields_num()-1);

    cout<<"[NOTA] compute gauss-cCurvature"<<endl;
    ccurvature = C_Curvature(GAUSS);
    time.start();
    ccurvature.compute(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] gauss-cCurvature computation: ");
    ccurvature.print_curvature_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
    Writer::write_mesh_curvature_VTK(out.str(),tree.get_mesh(),"gauss_cCurvature",tree.get_mesh().get_vertex(1).get_fields_num()-1);
}

void compute_curvature(PMRT_Tree& tree, cli_parameters &cli)
{
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);

    cout<<"[NOTA] compute concentrated-curvature"<<endl;
    Timer time;
    Concentrated_Curvature curvature = Concentrated_Curvature(CONCENTRATED);
    time.start();
    curvature.compute_borders(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] border computation: ");
    time.start();
    curvature.compute(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] concentrated-curvature computation: ");
    curvature.print_curvature_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
    Writer::write_mesh_curvature_VTK(out.str(),tree.get_mesh(),"concentrated",tree.get_mesh().get_vertex(1).get_fields_num()-1);

    cout<<"[NOTA] compute mean-cCurvature"<<endl;
    C_Curvature ccurvature = C_Curvature(MEAN,true);
    time.start();
    ccurvature.compute(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] mean-cCurvature computation: ");
    ccurvature.print_curvature_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
    Writer::write_mesh_curvature_VTK(out.str(),tree.get_mesh(),"mean_cCurvature_sumAll",tree.get_mesh().get_vertex(1).get_fields_num()-1);

    cout<<"[NOTA] compute gauss-cCurvature"<<endl;
    ccurvature = C_Curvature(GAUSS);
    time.start();
    ccurvature.compute(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] gauss-cCurvature computation: ");
    ccurvature.print_curvature_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
    Writer::write_mesh_curvature_VTK(out.str(),tree.get_mesh(),"gauss_cCurvature",tree.get_mesh().get_vertex(1).get_fields_num()-1);
}
