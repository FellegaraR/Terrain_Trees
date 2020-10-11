#ifndef GRAD_AWARE_CONTRACTION_H
#define GRAD_AWARE_CONTRACTION_H

#include "contraction_simplifier.h"
#include "morse/forman_gradient.h"
#include "morse/forman_gradient_topological_relations.h"

class Gradient_Aware_Simplifier: Contraction_Simplifier{


public:
    Gradient_Aware_Simplifier(){}
    void gradient_aware_simplify(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli,Forman_Gradient &gradient);
protected:
    void simplify_compute(Node_V &n,  Mesh &mesh, LRU_Cache<int, leaf_VT> &cache,Spatial_Subdivision &division,  contraction_parameters &params,PRT_Tree &tree,Forman_Gradient &gradient);
    void simplify_leaf(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, PRT_Tree& tree,Forman_Gradient &gradient); 
    void simplify_leaf_QEM(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, PRT_Tree& tree,Forman_Gradient &gradient);
    void contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1, Node_V &outer_v_block, edge_queue &edges,
                              Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params,Forman_Gradient &gradient);
    // void get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1, Node_V *& outer_v_block, 
    //                                Node_V &n, Mesh &mesh, leaf_VT &vts, LRU_Cache<int,leaf_VT> &cache, contraction_parameters &params, PRT_Tree &tree);
    bool valid_gradient_configuration(int v1,int v2,VT &vt1, VT &vt2,ET& et, bool v1_is_border, bool v2_is_border, Forman_Gradient &gradient, Mesh &mesh);
    void update(const ivect &e, VT& vt, VT& difference, Node_V &n, Node_V &v_block, edge_queue &edges,
                                          Mesh &mesh, contraction_parameters &params);
};

#endif // GRAD_AWARE_CONTRACTION_H