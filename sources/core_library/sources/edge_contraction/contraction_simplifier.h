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
    static void simplify_leaf(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params);

    /// skips an edge that has one of the two extreme deleted
    inline static bool skip_edge(const ivect &e, Mesh &mesh) /// the procedure checks if one of the two extreme is already deleted
    {
        return (mesh.is_vertex_removed(e[0]) || mesh.is_vertex_removed(e[1]));
    }

    ///
    static void contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1, leaf_VT &vtops, Node_V &outer_v_block, edge_queue &edges,
                              Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params);
    /// initializes the VTop and ETop of the edge and its two vertices
    /// plus saves the second leaf block if we are processing a cross-edge
    template<class T> static void get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1, Node_V *& outer_v_block,
                                   Node_V &n, Mesh &mesh, leaf_VT &vtops, LRU_Cache<int,leaf_VT> &cache, contraction_parameters &params, T &tree);
    /// the VTop is always without removed top-simplices
    template<class T> static VT* get_VT(int v_id, Node_V &n, Mesh &mesh, leaf_VT &vts, LRU_Cache<int,leaf_VT> &cache,
                        T &tree, Node_V *& v_block, contraction_parameters &params);
     // Find two adjacent triangles of edge e.                   
    static void get_ET(ivect &e, ET &et, Node_V &n, Mesh &mesh, leaf_VT &vts);
    ///

    /// the procedure removes from the VT relation the deleted top d-simplices
    static void clean_coboundary(VT &cob, Mesh &mesh);


    /// the procedure updates
    /// (1) the VTop relation of the surviving vertex
    /// (2) the top simplices for which we change the v2 with v1
    /// (3) checks and updates the top list in the target leaf block n (if a top simplex is now incident in it)
    /// (4) and add the new edges to the edge queue
    static void update(const ivect &e, VT& vt, VT& difference, Node_V &n, Node_V &v_block, edge_queue &edges,
                                          Mesh &mesh, contraction_parameters &params);
    static void remove_from_mesh(int to_delete_v, ET &et, Mesh &mesh, contraction_parameters &params);    
};

#endif // CONTRACTION_SIMPLIFIER_H
