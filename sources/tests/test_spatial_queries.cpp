#include "utilities/utility_functions.h"

using namespace utility_functions;

template<class T> void load_tree(T& tree, cli_parameters &cli);

int main(int , char** )
{
    cerr<<"[OBJECTIVE] this unit-test generates a quadtrees based on the PR-T tree criterion. "
        <<"Then it computes the critical points, 10 point locations, 10 neighbor finding searches and 10 box queries."<<endl;
    cerr<<"[TEST] INPUT MESH eggs64z4.tri"<<endl;
    cli_parameters cli;
    cli.mesh_path = "../data/eggs64z4.tri";
	cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = 10;
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PR-T tree"<<endl;
    load_tree(ptree,cli);

    Timer time;
    cout<<"[TEST] Critical Points extraction"<<endl;
    Critical_Points_Extractor cpe = Critical_Points_Extractor();
    time.start();
    cpe.compute_critical_points(ptree.get_root(),ptree.get_mesh(),ptree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] critical points extraction: ");
    cpe.print_stats();

    cout<<"[TEST] Random Points generation"<<endl;
    set<Point> points;
    Input_Generator::generate_random_point_inputs(ptree.get_mesh().get_domain(),10,points);

    Statistics stats;
    Spatial_Queries sq;
    cout<<"[TEST] Executing Point Locations"<<endl;
    sq.exec_point_locations(ptree,points,stats);
    cout<<"[TEST] Executing Nearest Neighbor Searches -- MAXIMA"<<endl;
    sq.exec_incremental_nearest_neighbor_queries(ptree,points,Point_Type::MAXIMUM,cpe.get_critical_points());

    stats = Statistics();
    points.clear();

    cout<<"[TEST] Random Boxes generation (5% of diagonal length)"<<endl;
    set<Box> boxes;
    Input_Generator::generate_random_box_inputs(ptree.get_mesh().get_domain(),0.05,10,boxes);

    cout<<"[TEST] Executing Box Queries"<<endl;
    sq.exec_box_queries(ptree,boxes,stats);


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
    
    time.start();
    Reindexer reindexer = Reindexer();
    reindexer.reindex_tree_and_mesh(tree,false,cli.original_vertex_indices,
                                    false,cli.original_triangle_indices);
    time.stop();
    time.print_elapsed_time("[TIME] Index and Mesh Reindexing ");    
}
