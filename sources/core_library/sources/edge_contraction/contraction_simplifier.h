#ifndef CONTRACTION_SIMPLIFIER_H
#define CONTRACTION_SIMPLIFIER_H


#include "terrain_trees/node_v.h"
#include "terrain_trees/node_t.h"
#include "utilities/lru_cache.h"
#include "terrain_trees/prt_tree.h"
//#include "topological_ds/links_aux_structures.h"  SHOULD BE REPLACED WITH WHICH HEADER?
#include "statistics/statistics.h"
#include "utilities/container_utilities.h"
#include "terrain_trees/mesh_updater.h"
#include "utilities/usage.h"
#include "utilities/cli_parameters.h"
#include "utilities/string_management.h"
//#include "topological_queries/topological_query_extractor.h"
#include "queries/topological_queries.h"

#include "simplification_aux_structures.h"



/// Contraction_Simplifier class ///
class Contraction_Simplifier
{
public:
    Contraction_Simplifier() {}
    /// this procedure simplify the simplicial complex without considering any weight for the edges
    template<class T> void simplify(T &tree, Mesh &mesh, cli_parameters &cli);
    template<class T> static void update_mesh_and_tree(T &tree, Mesh &mesh,contraction_parameters &params);

protected:
    template<class T> static void simplify_compute(Node_V &n,  Mesh &mesh, LRU_Cache<int, leaf_VT> &cache,Spatial_Subdivision &division,  contraction_parameters &params, T &tree);
    template<class T>   static void simplify_leaf(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, T& tree);

    /// skips an edge that has one of the two extreme deleted
    inline static bool skip_edge(const ivect &e, Mesh &mesh) /// the procedure checks if one of the two extreme is already deleted
    {
        return (mesh.is_vertex_removed(e[0]) || mesh.is_vertex_removed(e[1]));
    }

    ///
    static void contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1, Node_V &outer_v_block, edge_queue &edges,
                              Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params);
    /// initializes the VTop and ETop of the edge and its two vertices
    /// plus saves the second leaf block if we are processing a cross-edge
    template<class T> static void get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1, Node_V *& outer_v_block,
                                   Node_V &n, Mesh &mesh, leaf_VT &vts, LRU_Cache<int,leaf_VT> &cache, contraction_parameters &params, T &tree);
    /// the VTop is always without removed top-simplices
    template<class T> static VT* get_VT(int v_id, Node_V &n, Mesh &mesh, leaf_VT &vts, LRU_Cache<int,leaf_VT> &cache,
                        T &tree, Node_V *& v_block, contraction_parameters &params);
     // Find two adjacent triangles of edge e.                   
    static void get_ET(ivect &e, ET &et, Node_V &n, Mesh &mesh, leaf_VT &vts);
    ///

    /// the procedure removes from the VT relation the deleted top d-simplices
    static void clean_coboundary(VT &cob, Mesh &mesh);

    static void find_candidate_edges(Node_V &n, Mesh &mesh, leaf_VT &vts,edge_queue& edges, contraction_parameters &params);
    /// the procedure updatess
    /// (1) the VTop relation of the surviving vertex
    /// (2) the top simplices for which we change the v2 with v1
    /// (3) checks and updates the top list in the target leaf block n (if a top simplex is now incident in it)
    /// (4) and add the new edges to the edge queue
    static void update(const ivect &e, VT& vt, VT& difference, Node_V &n, Node_V &v_block, edge_queue &edges,
                                          Mesh &mesh, contraction_parameters &params);
    static void remove_from_mesh(int to_delete_v, ET &et, Mesh &mesh, contraction_parameters &params);   
    static bool link_condition(int v0, int v1, VT &vt0, VT &vt1, Mesh &mesh); 
};


template<class T> void Contraction_Simplifier::simplify(T &tree, Mesh &mesh, cli_parameters &cli)
{   

    cerr<<"==Homology preserving simplification - weak-link condition=="<<endl;

    cerr<<"[NOTICED] Cache size: "<<cli.cache_size<<endl;
    LRU_Cache<int,leaf_VT> cache(cli.cache_size); // the key is v_start while the value are the VTop relations
    contraction_parameters params;
    params.set_maximum_length(cli.maximum_length);

    Timer time;
    int simplification_round;
    int round = 1;

    time.start();
    while(1)
    {
        simplification_round = params.get_contracted_edges_num();  //checked edges
        /// HERE YOU NEED TO DEFINE A PROCEDURE FOR SIMPLIFY THE TIN BY USING THE SPATIAL INDEX
        this->simplify_compute(tree.get_root(),mesh,cache,tree.get_subdivision(),params,tree);
        // PARTIAL SIMPLIFICATION STATS
        cerr<<"=== end-of-round "<<round<<") --> contracted edges: ";
        cerr<<params.get_contracted_edges_num()-simplification_round<<endl;
        round++;

        if(cli.debug_mode)
        {
            time.stop();
            time.print_elapsed_time("   [TIME] executing a simplification round: ");
            time.start();
          //  cerr << "   [RAM] peak for executing a simplification round: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
        }

        if(simplification_round == params.get_contracted_edges_num())
            break;

        cache.reset();
    }
    time.stop();
    if(!cli.debug_mode)
        time.print_elapsed_time("[TIME] Edge contraction simplification: ");
    //else
        //params.print_simplification_partial_timings();
    //params.print_simplification_counters();

  //  cerr << "[RAM peak] for contracting a simplicial complex: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

    /// finally we have to update/compress the mesh and the tree
    Contraction_Simplifier::update_mesh_and_tree(tree,mesh,params);
}

  template<class T>  void Contraction_Simplifier::simplify_compute(Node_V &n,  Mesh &mesh, LRU_Cache<int, leaf_VT> &cache,Spatial_Subdivision &division,  contraction_parameters &params, T &tree)
{
      if (n.is_leaf())
    {
        simplify_leaf(n,mesh,cache,params,tree);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                simplify_compute(*n.get_son(i),mesh,cache,division,params,tree);
            }
        }
    }

}



template<class T> void Contraction_Simplifier::simplify_leaf(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, T& tree){

if(!n.indexes_vertices())
     return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();
    itype v_range = v_end - v_start;

leaf_VT local_vts(v_range,VT());
n.get_VT(local_vts,mesh);
// Create a priority quue of candidate edges
edge_queue edges;
find_candidate_edges(n,mesh,local_vts,edges,params);
    while(!edges.empty())
    {
        Geom_Edge* current = edges.top();
        edges.pop();
        ET et;
        VT *vt0,*vt1;
        Node_V *outer_v_block;
        //leaf_VT vts;
         ivect e=current->edge;
        get_edge_relations(e,et,vt0,vt1,outer_v_block,n,mesh,local_vts,cache,params,tree);
        if(link_condition(e[0],e[1],*vt0,*vt1,mesh)){
        contract_edge(e,et,*vt0,*vt1,*outer_v_block,edges,n,mesh,cache,params);
    
        }
    }


}

template<class T> void Contraction_Simplifier::get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1, Node_V *&outer_v_block,
                                                Node_V &n, Mesh &mesh, leaf_VT &vts, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params,T &tree)
{


    outer_v_block = NULL;
    /// inverted order as I only need the block indexing v1
    vt1 = Contraction_Simplifier::get_VT(e[1],n,mesh,vts,cache,tree,outer_v_block,params);
    vt0 = Contraction_Simplifier::get_VT(e[0],n,mesh,vts,cache,tree,outer_v_block,params);



    Contraction_Simplifier::get_ET(e,et,n,mesh,vts);


}

template<class T>  VT* Contraction_Simplifier::get_VT(int v_id, Node_V &n, Mesh &mesh, leaf_VT &vts, LRU_Cache<int,leaf_VT> &cache,
                        T &tree, Node_V *& v_block, contraction_parameters &params)
{
    int local_index;
    bool debug=false;
    if(n.indexes_vertex(v_id))
    {
        if(debug)
            cout<<"[get_VT] "<<v_id<<" -> LOCAL VERTEX "<<n<<endl;

        local_index = v_id - n.get_v_start();
        VT* vt = &(vts[local_index]);
        Contraction_Simplifier::clean_coboundary(*vt,mesh);
        v_block = &n;
        return vt;
    }
    else
    {
        if(debug)  // if v is external
            cout<<"[get_VT] "<<v_id<<" -> EXTERNAL VERTEX "<<n<<endl;

        tree.get_leaf_indexing_vertex(tree.get_root(),v_id,v_block);
        local_index = v_id - v_block->get_v_start();

        LRU_Cache<int,leaf_VT>::mapIt it_c = cache.find(v_block->get_v_start()); //First check in the cache
        if(it_c == cache.end())   //if not in the cache
        {
            if(debug)
                cout<<"    -> LEAF BLOCK OUTSIDE CACHE - REGEN "<<*v_block<<endl;

            leaf_VT lVT;
            v_block->get_VT(lVT,mesh);
            it_c = cache.insert(v_block->get_v_start(),lVT);
        }
        else
        {
            if(debug)
            {
                cout<<"    -> LEAF BLOCK IN CACHE - CLEAN "<<*v_block<<endl;
            }

            // if(debug/* || v_id ==2355*/)
            //     cout<<"num_elem_in_vt: "<<get_num_elements_in_container_of_containers((it_c->second)[local_index])<<endl;
//                print_container_of_containers_content("VTop(2355) ",(it_c->second)[local_index]);

            Contraction_Simplifier::clean_coboundary((it_c->second)[local_index],mesh);

            // if(debug/* || v_id ==2355*/)
            //     cout<<"num_elem_in_vt: "<<get_num_elements_in_container_of_containers((it_c->second)[local_index])<<endl;
//                print_container_of_containers_content("CleanedVTop(2355) ",(it_c->second)[local_index]);
        }

        return &(it_c->second)[local_index];
    }
}


template<class T> void Contraction_Simplifier::update_mesh_and_tree(T &tree, Mesh &mesh, contraction_parameters &params)
{
    Timer time;

    ///  UPDATE OF MESH AND TREE
    ivect new_v_positions;
    ivect new_t_positions;
    ivect surviving_vertices;

    time.start();
//    cerr<<"[TREE] compact vertices lists"<<endl;
    tree.compact_vertices_lists(tree.get_root(),mesh,surviving_vertices);
    time.stop();
    time.print_elapsed_time("[TIME] Compact tree vertices lists: ");

//    print_container_content("surviving vertices: ",surviving_vertices);
//    mesh.print_mesh(cout);
//    int a; cin>>a;

    time.start();
//    cerr<<"[MESH] compact"<<endl;
    Mesh_Updater mu;
    mu.clean_vertices_array(mesh,new_v_positions,surviving_vertices);
    /// NEW: the update_and_compact procedure check internally if we have removed all the top d-simplices
    boost::dynamic_bitset<> all_deleted = mu.update_and_clean_triangles_arrays(mesh,new_v_positions,new_t_positions,params.get_counter());
    time.stop();
    time.print_elapsed_time("[TIME] Compact and update mesh: ");

    // cerr<<"[STAT] mesh "<<endl;
    // cerr<<"------------"<<endl;
    // mesh.print_mesh_stats(cerr);
    // cerr<<"------------"<<endl;

    time.start();
//    cerr<<"[TREE] update indices in the tree"<<endl;
    tree.update_tree(tree.get_root(),new_v_positions,new_t_positions,all_deleted);
    time.stop();
    time.print_elapsed_time("[TIME] Update tree (top-simplices): ");

//    Reindexer r;
//    r.reorganize_index_and_mesh(tree,mesh,false);

    //cerr << "[RAM peak] for updating the mesh and the tree: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
}


#endif // CONTRACTION_SIMPLIFIER_H
