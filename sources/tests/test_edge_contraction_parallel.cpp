#include "utilities/utility_functions.h"

using namespace utility_functions;

void load_tree(PRT_Tree &tree, cli_parameters &cli);

int main(int argc, char **argv)
{
  cli_parameters cli;
  cli.mesh_path = argv[1];

  cerr << "[OBJECTIVE] this unit-test generates a PR-quadtree on the input TIN dataset "
       << "then, it simplifies the triangle mesh with an edge contraction operator following a length criteria." << endl;

  cerr << "[NOTA] all the generated files are saved in the 'data' folder" << endl;

  cli.division_type = QUAD;
  cli.crit_type = "pr";
  cli.v_per_leaf = atoi(argv[2]);
      if(atof(argv[4])==-1){
        cli.contract_all_edges = true;
    }
  
  cli.maximum_limit = atof(argv[4]);
 if(atoi(argv[5])==-1){
   cli.num_of_threads = omp_get_max_threads();
 }
 else{
  cli.num_of_threads =atoi(argv[5]);
 }
  cli.debug_mode = false;
  PRT_Tree ptree = PRT_Tree(cli.v_per_leaf, cli.division_type);
  cerr << "[GENERATION] PR-T tree" << endl;
  if (strcmp(argv[3], "-q") == 0)
    cli.QEM_based = true;
  else
    cli.QEM_based = false;
  load_tree(ptree, cli);

  return (EXIT_SUCCESS);
}

void load_tree(PRT_Tree &tree, cli_parameters &cli)
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

    stringstream out;

    out << base.str() << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type;
    out << "_v_" << cli.v_per_leaf << "_.tree";
    cli.tree_path=out.str();
  if (!Reader::read_tree(tree, tree.get_root(), cli.tree_path))
  {
 
    cerr << "[NOTICE] Cannot find existed .tree file." << endl;
    cerr << "[GENERATION] tree from triangle mesh" << endl;
    time.start();
    tree.build_tree();
    time.stop();
    time.print_elapsed_time(tree_info.str());
    Writer::write_tree(out.str(), tree.get_root(), tree.get_subdivision());
  }
  else
    cout << "[NOTICE] Found corresponding .tree file. Loaded tree from file successfully" << endl;

  cout << "number of triangles: " << tree.get_mesh().get_triangles_num() << endl;
  cout << "number of vertices: " << tree.get_mesh().get_vertices_num() << endl;
 
  stringstream out2;
  out2 << base.str();
  out2 << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type << "_v_" << cli.v_per_leaf << "_tree.vtk";

  string output_name = base.str() + "_simplified_v_" + to_string(cli.v_per_leaf);
  if (cli.QEM_based == true)
  {
    output_name = output_name + "_qem";
  }
  else
  {
    output_name = output_name + "_length";
  }
  Writer::write_tree_VTK(out2.str(), tree.get_root(), tree.get_subdivision(), tree.get_mesh());
  // Writer::write_mesh_VTK(base.str(),tree.get_mesh());
  
  time.start();
  Reindexer reindexer = Reindexer();
  reindexer.reindex_tree_and_mesh(tree, false, cli.original_vertex_indices,
                                  false, cli.original_triangle_indices);
  time.stop();
  time.print_elapsed_time("[TIME] Index and Mesh Reindexing ");

  tree.init_leaves_list(tree.get_root()); // NOTA: enables the parallel OpenMP-based processing of the leaf block

  Statistics stats;
  stats.get_index_statistics(tree, cli.reindex);

  cout << "[NOTA]Border checking" << endl;
  time.start();
  Border_Checker border_checker = Border_Checker();
  border_checker.compute_borders(tree.get_root(), tree.get_mesh().get_domain(), 0, tree.get_mesh(), tree.get_subdivision());

  time.stop();
  time.print_elapsed_time("[TIME] Border Checking ");
  time.start();
  Contraction_Simplifier simplifier;
  simplifier.preprocess(tree, tree.get_mesh(), cli);
  time.stop();
  time.print_elapsed_time("[TIME] Preporcessing");
  time.start();
  simplifier.simplify_parallel(tree, tree.get_mesh(), cli);
  time.stop();
  time.print_elapsed_time("[TIME] Complete edge contraction process ");
  cout << "number of remaining triangles: " << tree.get_mesh().get_triangles_num() << endl;

  cout << output_name << endl;
  Writer::write_mesh_VTK(output_name, tree.get_mesh());
  Writer::write_mesh(output_name, "parallel", tree.get_mesh(), false);
}
