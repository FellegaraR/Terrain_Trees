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

#define _v0 0
#define _v1 1


/// Contraction_Simplifier class ///
class Contraction_Simplifier
{
public:
    Contraction_Simplifier() {}
    /// this procedure simplify the simplicial complex without considering any weight for the edges
    void simplify(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli);
    void update_mesh_and_tree(PRT_Tree &tree, Mesh &mesh,contraction_parameters &params);
    

protected:
    
   void simplify_compute(Node_V &n,  Mesh &mesh, LRU_Cache<int, leaf_VT> &cache,Spatial_Subdivision &division,  contraction_parameters &params,PRT_Tree &tree);
    void simplify_leaf(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, PRT_Tree& tree);
    ///edge contraction based on QEM criterion 
    void simplify_leaf_QEM(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, PRT_Tree& tree);
    /// skips an edge that has one of the two extreme deleted
    inline bool skip_edge(const ivect &e, Mesh &mesh) /// the procedure checks if one of the two extreme is already deleted
    {
        return (mesh.is_vertex_removed(e[0]) || mesh.is_vertex_removed(e[1]));
    }

    ///
     void contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1, Node_V &outer_v_block, edge_queue &edges,
                              Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params);
    /// initializes the VTop and ETop of the edge and its two vertices
    /// plus saves the second leaf block if we are processing a cross-edge
    void get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1, Node_V *& outer_v_block,
                                   Node_V &n, Mesh &mesh, leaf_VT &vts, LRU_Cache<int,leaf_VT> &cache, contraction_parameters &params, PRT_Tree &tree);
    /// the VTop is always without removed top-simplices
     VT* get_VT(int v_id, Node_V &n, Mesh &mesh, leaf_VT &vts, LRU_Cache<int,leaf_VT> &cache,
                        PRT_Tree &tree, Node_V *& v_block, contraction_parameters &params);
     // Find two adjacent triangles of edge e.                   
    void get_ET(ivect &e, ET &et, Node_V &n, Mesh &mesh, leaf_VT &vts);
    ///

    /// the procedure removes from the VT relation the deleted top d-simplices
    void clean_coboundary(VT &cob, Mesh &mesh);

    void find_candidate_edges_QEM(Node_V &n, Mesh &mesh, leaf_VT &vts,edge_queue& edges, contraction_parameters &params);
    void find_candidate_edges(Node_V &n, Mesh &mesh, leaf_VT &vts,edge_queue& edges, contraction_parameters &params);
    /// the procedure updatess
    /// (1) the VTop relation of the surviving vertex
    /// (2) the top simplices for which we change the v2 with v1
    /// (3) checks and updates the top list in the target leaf block n (if a top simplex is now incident in it)
    /// (4) and add the new edges to the edge queue
    void update(const ivect &e, VT& vt, VT& difference, Node_V &n, Node_V &v_block, edge_queue &edges,
                                          Mesh &mesh, contraction_parameters &params);
    void remove_from_mesh(int to_delete_v, ET &et, Mesh &mesh, contraction_parameters &params);   
    bool link_condition(int v0, int v1, VT &vt0, VT &vt1, Mesh &mesh); 
    void update_new(const ivect &e, VT& vt, VT& difference, Node_V &n, Node_V &v_block, edge_queue &edges,
                                          Mesh &mesh, contraction_parameters &params, int new_vertex_pos);
    void compute_initial_QEM(Mesh &mesh, vector<dvect >& planes);
    void compute_triangle_plane( Mesh &mesh,vector<dvect>& trPl);
   // new_vertex is the pos of the remaining vertex in the edge to be contracted. 0 refers to e[0], 1 refers to e[1]
    double compute_error(int v1, int v2, Mesh &mesh,  int& new_vertex_pos);
    inline double vertex_error(Matrix q, double x, double y, double z)
{
    return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x + q[5]*y*y
        + 2*q[6]*y*z + 2*q[7]*y + q[10]*z*z + 2*q[11]*z + q[15];
}

      vector<Matrix> initialQuadric;
      vector<dvect> trianglePlane;
      map<vector<int>,double> updated_edges;
};



#endif // CONTRACTION_SIMPLIFIER_H
