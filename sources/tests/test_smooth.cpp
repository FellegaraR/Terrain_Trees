#include "utilities/utility_functions.h"

using namespace utility_functions;

void load_tree(PRT_Tree& tree, cli_parameters &cli);
void smooth_mesh(PRT_Tree& tree, cli_parameters &cli, double lambda, int iter_num);


int main(int argc , char** argv)
{
	cli_parameters cli;
    cli.mesh_path = argv[1];
    
    cerr<<"[OBJECTIVE] this unit-test generates a quadtrees based on the PR-T tree criterion. "
	    <<"then, it smoothes all the vertices in the mesh. "<<endl;
		
    cerr<<"[NOTA] all the generated files are saved in the 'data' folder"<<endl;
    
    cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = atoi(argv[2]);    
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PR-T tree"<<endl;
    load_tree(ptree,cli);
    smooth_mesh(ptree,cli,atof(argv[3]),atoi(argv[4]));

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
    stringstream tree_info;
    tree_info << base_info.str() << "[TIME] Building ";
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


    stringstream out2;
    out2 << base.str();
    out2 << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type << "_v_" << cli.v_per_leaf << "_tree.vtk";
    


    time.start();
    Reindexer reindexer = Reindexer();
    reindexer.reindex_tree_and_mesh(tree,false,cli.original_vertex_indices,
                                    false,cli.original_triangle_indices);
    time.stop();
    time.print_elapsed_time("[TIME] Index and Mesh Reindexing ");
         cerr << "[MEMORY] peak for Index and Mesh Reindexing: " <<
        to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
}  

void smooth_mesh(PRT_Tree& tree,  cli_parameters &cli, double lambda, int iter_num)
{
    stringstream output_name;
    output_name << get_path_without_file_extension(cli.mesh_path);

    tree.init_leaves_list(tree.get_root()); 
    cout<<"[NOTA]Border checking"<<endl;
    Border_Checker border_checker=Border_Checker();
    border_checker.compute_borders(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    cout<<"[NOTA] compute roughness"<<endl;
    Timer time;
    Smooth smoother = Smooth();
    time.start();
    smoother.smooth_parallel(tree,tree.get_mesh(),lambda,iter_num,cli);
    time.stop();
    time.print_elapsed_time("[TIME] smooothing computation: ");
    cerr << "[MEMORY] peak for smooothing: " <<
    to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
 
    Writer::write_mesh(output_name.str(),"smoothed",tree.get_mesh(),false); 
   }

