#include "utilities/utility_functions.h"
//#include "utility_functions.h"

using namespace utility_functions;

template<class T> void load_tree(T& tree, cli_parameters &cli);
template<class T> void exec_queries(T& tree, cli_parameters &cli);
template<class T> void compute_curvature(T& tree, cli_parameters &cli);
void compute_curvature(PMRT_Tree& tree, cli_parameters &cli);
template<class T> void compute_terrain_features(T& tree, cli_parameters &cli);
void compute_terrain_features(PMRT_Tree& tree, cli_parameters &cli);
void generate_query_inputs(cli_parameters &cli);

int main(int argc, char** argv)
{
    if(argc == 1)
    {
        print_help();
        return 0;
    }

    cli_parameters cli = cli_parameters();

    if (read_arguments(argc, argv,cli) == -1)
    {
        print_usage();
        return (EXIT_FAILURE);
    }

    if (cli.isTreeFile)
    {
        //nel caso di lettura da file recupero le info dal nome del file
        setParameters(cli);
    }

    //controllo che tutto sia inizializzato correttamente
    if (!cli.is_getInput)
    {
        if (!checkParameters(cli))
            return (EXIT_FAILURE);
    }

    if (cli.crit_type == "pr")
    {
        PRT_Tree tree = PRT_Tree(cli.v_per_leaf,cli.division_type);
        load_tree(tree,cli);
        exec_queries(tree,cli);
        if(cli.reindex)
        {
            compute_curvature(tree,cli);            
            compute_terrain_features(tree,cli);
        }
    }
    else if (cli.crit_type == "pm")
    {
        PMT_Tree tree = PMT_Tree(cli.v_per_leaf, cli.t_per_leaf,cli.division_type);
        load_tree(tree,cli);
        exec_queries(tree,cli);
        if(cli.reindex)
        {
            compute_curvature(tree,cli);            
            compute_terrain_features(tree,cli);
        }
    }
    else if (cli.crit_type == "pmr")
    {
        PMRT_Tree tree = PMRT_Tree(cli.t_per_leaf,cli.division_type);
        load_tree(tree,cli);
        exec_queries(tree,cli);
        if(cli.reindex)
        {
            compute_curvature(tree,cli);
            compute_terrain_features(tree,cli);            
        }
    }
    else if (cli.is_getInput)
    {
        generate_query_inputs(cli);
    }

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
    else if(cli.query_type == FILTER)
    {
        ///read a point cloud
        if(!Reader::read_vertices(tree.get_mesh(),cli.mesh_path))
        {
            cerr << "[ERROR] Loading point cloud file. Execution Stopped." << endl;
            return;
        }
    }
    /// tri or off format triangle mesh
    else if (!Reader::read_mesh(tree.get_mesh(), cli.mesh_path))
    {
        cout << "[ERROR] Loading mesh file. Execution Stopped." << endl;
        return;
    }

    cerr << "[MEMORY] peak for Indexing the terrain: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    stringstream base_info;
    base_info << cli.v_per_leaf << " " << cli.t_per_leaf << " " << cli.crit_type << " ";
    stringstream base;
    base << get_path_without_file_extension(cli.mesh_path);

    stringstream tree_info;
    tree_info << base_info.str();

    /// if we want to filter a point clouds we force the construction of the tree
    if (cli.isTreeFile && cli.query_type != FILTER)
    {
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
        time.start();
        if(get_file_extension(cli.mesh_path) == "soup")
        {
            cerr<<"[GENERATION] tree from soup of triangles"<<endl;
            tree.build_tree(soup);
            soup.clear();
            Writer::write_mesh(base.str(),"from_soup",tree.get_mesh(),false);
        }
        else if(cli.query_type == FILTER)
        {
            /// here we generate a PR tree on the points
            /// we filter them and finally we output two files containing
            /// the clouds in a SpatialHadoop can understand
            /// and a multifield version of the cloud for future usage
            cerr<<"[GENERATION] tree from points cloud"<<endl;
            vertex_multifield multifield;
            tree.build_tree_from_cloud(multifield);
            time.stop();
            time.print_elapsed_time(tree_info.str());
            cerr<<"[STATS] vertices: "<<tree.get_mesh().get_vertices_num()<<" vs. multifield_vertices: "<<multifield.size()<<endl;
            cerr<<"[OUTPUT] writing points cloud files"<<endl;
            Writer::write_filtered_points_cloud(base.str(),tree.get_mesh());
            Writer::write_multifield_points_cloud(base.str(),multifield,tree.get_mesh());
            return; /// after generated these files we have done..
        }
        else
        {
            cerr<<"[GENERATION] tree from triangle mesh"<<endl;
            tree.build_tree();
        }
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

        if(cli.app_debug == OUTPUT)
        {
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
        }
    }

    cerr << "[MEMORY] peak for encoding the Terrain tree: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    if(cli.query_type == MORSE_ANALYSIS)
    {
        cli.rever_to_original = Reader::read_noised_field_value(tree.get_mesh(),base.str(),cli.original_vertex_fields);
    }

    if(cli.reindex)
    {
        cerr<<"[REINDEXING] tree and triangle mesh"<<endl;

        if(cli.query_type == MORSE_ANALYSIS || cli.query_type == LOCAL_MORSE_SIMPLIFICATION || cli.query_type == GLOBAL_MORSE_SIMPLIFICATION)
            cli.original_vertex_indices.assign(tree.get_mesh().get_vertices_num(),-1);
        if(cli.query_type == MORSE_ANALYSIS && cli.app_debug == OUTPUT)
            cli.original_triangle_indices.assign(tree.get_mesh().get_triangles_num(),-1);

        time.start();
        Reindexer reindexer = Reindexer();
        reindexer.reindex_tree_and_mesh(tree,(cli.query_type == MORSE_ANALYSIS),cli.original_vertex_indices,
                                        (cli.query_type == MORSE_ANALYSIS && cli.app_debug == OUTPUT),cli.original_triangle_indices);
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

template<class T> void exec_queries(T& tree, cli_parameters &cli)
{
    if (cli.query_type != NULL_QUERY)
    {
        stringstream base_info;
        base_info << cli.v_per_leaf << " " << cli.t_per_leaf << " " << cli.crit_type << " ";
        Statistics stats;
        Spatial_Queries sq;
        Topological_Queries tq;

        cerr<<base_info.str()<<endl;
        if (cli.query_type == POINT)
            sq.exec_point_locations(tree,cli.query_path,stats);
        else if(cli.query_type == BOX)
            sq.exec_box_queries(tree,cli.query_path,stats);
        else if(cli.query_type == WINDVT)
            tq.windowed_VT(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_subdivision(),cli.query_path,cli.reindex);
        else if(cli.query_type == WINDTT)
            tq.windowed_TT(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_subdivision(),cli.query_path);
        else if(cli.query_type == BATCH)
        {
            tq.batched_VT(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_subdivision(),cli.reindex);
            cerr << "[MEMORY] peak for extracting batched VT relations: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
            tq.batched_TT(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
            cerr << "[MEMORY] peak for extracting batched TT relations: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
        }
    }
}

template<class T> void compute_curvature(T& tree, cli_parameters &cli)
{    
    if(cli.query_type != CONCENTRATED_CURVATURE &&
            cli.query_type != MEAN_CCURVATURE &&
            cli.query_type != GAUSS_CCURVATURE)
        return;

    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);
    Timer time;

    if(cli.query_type == CONCENTRATED_CURVATURE)
    {
        cout<<"[NOTA] compute concentrated-curvature"<<endl;
        Concentrated_Curvature curvature = Concentrated_Curvature(CONCENTRATED);
        time.start();
        curvature.compute_borders(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] border computation: ");
        cerr << "[MEMORY] peak for computing the terrain borders: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
        time.start();
        curvature.compute(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] concentrated-curvature computation: ");
        curvature.print_curvature_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
        cerr << "[MEMORY] peak for computing the Concentrated Curvature: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
        if(cli.app_debug==OUTPUT)
            Writer::write_mesh_curvature_VTK(out.str(),tree.get_mesh(),"concentrated",tree.get_mesh().get_vertex(1).get_fields_num()-1);
    }

    if(cli.query_type == MEAN_CCURVATURE)
    {
        cout<<"[NOTA] compute mean-cCurvature"<<endl;
        C_Curvature ccurvature = C_Curvature(MEAN,true/*,(cli.app_debug == TIME_VERBOSE)*/);
        time.start();
        ccurvature.compute(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] mean-cCurvature computation: ");
        ccurvature.print_curvature_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
        cerr << "[MEMORY] peak for computing the Mean CCurvature: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
        if(cli.app_debug==OUTPUT)
            Writer::write_mesh_curvature_VTK(out.str(),tree.get_mesh(),"mean_cCurvature_sumAll",tree.get_mesh().get_vertex(1).get_fields_num()-1);
//        if(cli.app_debug==TIME_VERBOSE)
//            ccurvature.print_debug_times();
    }

    if(cli.query_type == GAUSS_CCURVATURE)
    {
        cout<<"[NOTA] compute gauss-cCurvature"<<endl;
        C_Curvature ccurvature = C_Curvature(GAUSS,false/*,(cli.app_debug == TIME_VERBOSE)*/);
        time.start();
        ccurvature.compute(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] gauss-cCurvature computation: ");
        ccurvature.print_curvature_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
        cerr << "[MEMORY] peak for computing the Gauss CCurvature: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
        if(cli.app_debug==OUTPUT)
            Writer::write_mesh_curvature_VTK(out.str(),tree.get_mesh(),"gauss_cCurvature",tree.get_mesh().get_vertex(1).get_fields_num()-1);
//        if(cli.app_debug==TIME_VERBOSE)
//            ccurvature.print_debug_times();
    }
}

void compute_curvature(PMRT_Tree& tree, cli_parameters &cli)
{
    if(cli.query_type != CONCENTRATED_CURVATURE &&
            cli.query_type != MEAN_CCURVATURE &&
            cli.query_type != GAUSS_CCURVATURE)
        return;
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);
    Timer time;

    if(cli.query_type == CONCENTRATED_CURVATURE)
    {
        cout<<"[NOTA] compute concentrated-curvature"<<endl;
        Concentrated_Curvature curvature = Concentrated_Curvature(CONCENTRATED);        
        time.start();
        curvature.compute_borders(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] border computation: ");
        cerr << "[MEMORY] peak for computing the terrain borders: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
        time.start();
        curvature.compute(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] concentrated-curvature computation: ");
        curvature.print_curvature_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
        cerr << "[MEMORY] peak for computing the Concentrated Curvature: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
        if(cli.app_debug==OUTPUT)
            Writer::write_mesh_curvature_VTK(out.str(),tree.get_mesh(),"concentrated",tree.get_mesh().get_vertex(1).get_fields_num()-1);
    }

    if(cli.query_type == MEAN_CCURVATURE)
    {
        cout<<"[NOTA] compute mean-cCurvature"<<endl;
        C_Curvature ccurvature = C_Curvature(MEAN,true);
        time.start();
        ccurvature.compute(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] mean-cCurvature computation: ");
        ccurvature.print_curvature_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
        cerr << "[MEMORY] peak for computing the Mean CCurvature: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
        if(cli.app_debug==OUTPUT)
            Writer::write_mesh_curvature_VTK(out.str(),tree.get_mesh(),"mean_cCurvature_sumAll",tree.get_mesh().get_vertex(1).get_fields_num()-1);
    }

    if(cli.query_type == GAUSS_CCURVATURE)
    {
        cout<<"[NOTA] compute gauss-cCurvature"<<endl;
        C_Curvature ccurvature = C_Curvature(GAUSS);
        time.start();
        ccurvature.compute(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] gauss-cCurvature computation: ");
        ccurvature.print_curvature_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
        cerr << "[MEMORY] peak for computing the Mean CCurvature: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
        if(cli.app_debug==OUTPUT)
            Writer::write_mesh_curvature_VTK(out.str(),tree.get_mesh(),"gauss_cCurvature",tree.get_mesh().get_vertex(1).get_fields_num()-1);
    }
}

void generate_query_inputs(cli_parameters &cli)
{
    Mesh mesh;
    //Legge l'input
    if (!Reader::read_mesh(mesh, cli.mesh_path))
    {
        cout << "Error Loading mesh file. Execution Stopped." << endl;
        return;
    }

    if (cli.is_getInput)
    {
        stringstream out;
        out << get_path_without_file_extension(cli.mesh_path);

        if(cli.input_gen_type == "rand")
        {
            if(cli.query_type == POINT && cli.ratio == 0.0)
                Input_Generator::generate_random_point_inputs(mesh.get_domain(),cli.num_input_entries, out.str());
            else if(cli.query_type == BOX && cli.ratio > 0.0)
                Input_Generator::generate_random_box_inputs(mesh.get_domain(),cli.ratio,cli.num_input_entries, out.str());
        }
        else if(cli.input_gen_type == "near")
        {
            if(cli.query_type == POINT && cli.ratio == 0.0)
                Input_Generator::generate_near_point_inputs(mesh.get_domain(),cli.num_input_entries, mesh, out.str());
            else if(cli.query_type == BOX && cli.ratio > 0.0)
                Input_Generator::generate_near_box_inputs(mesh.get_domain(),cli.ratio,cli.num_input_entries, mesh, out.str());
        }
    }
}

template<class T> void compute_terrain_features(T& tree, cli_parameters &cli)
{
    Timer time;
    if(cli.query_type == SLOPES)
    {
        Slope_Extractor se;
        time.start();
        se.compute_triangles_slopes(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] triangle-slopes computation: ");
        se.print_slopes_stats(tree.get_mesh().get_triangles_num());
        se.reset_stats();
        cerr << "[MEMORY] peak for computing the triangle-slopes: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

        time.start();
        se.compute_edges_slopes(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] edge-slopes computation: ");
        se.print_slopes_stats();
        se.reset_stats();
        cerr << "[MEMORY] peak for computing the edge-slopes: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
    }

    if(cli.query_type == CRITICAL_POINTS)
    {
        Critical_Points_Extractor cpe = Critical_Points_Extractor();
        time.start();
        cpe.compute_critical_points(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] critical points extraction: ");
        cpe.print_stats();
        cerr << "[MEMORY] peak for computing the critical points: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
    }
}

void compute_terrain_features(PMRT_Tree& tree, cli_parameters &cli)
{
    Timer time;
    if(cli.query_type == SLOPES)
    {
        Slope_Extractor se;
        time.start();
        se.compute_triangles_slopes(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] triangle-slopes computation: ");
        se.print_slopes_stats(tree.get_mesh().get_triangles_num());
        se.reset_stats();
        cerr << "[MEMORY] peak for computing the triangle-slopes: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

        time.start();
        se.compute_edges_slopes(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] edge-slopes computation: ");
        se.print_slopes_stats();
        se.reset_stats();
        cerr << "[MEMORY] peak for computing the edge-slopes: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
    }

    if(cli.query_type == CRITICAL_POINTS)
    {
        Critical_Points_Extractor cpe = Critical_Points_Extractor();
        time.start();
        cpe.compute_critical_points(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] critical points extraction: ");
        cpe.print_stats();
        cerr << "[MEMORY] peak for computing the critical points: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
    }
}
