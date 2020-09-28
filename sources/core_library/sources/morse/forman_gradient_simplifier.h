/*
    This file is part of the Terrain Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The Terrain Trees library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Terrain Trees library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Terrain Trees library.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FORMAN_GRADIENT_SIMPLIFICATION_H
#define FORMAN_GRADIENT_SIMPLIFICATION_H

#include "forman_gradient_features_extractor.h"
#include "forman_gradient_topological_relations.h"
#include "utilities/sorting.h"

typedef priority_queue<Topo_Sempl, vector<Topo_Sempl>, sort_arcs_topo> priority_arcs_queue;

class Forman_Gradient_Simplifier : public Forman_Gradient_Features_Extractor
{
public:
    Forman_Gradient_Simplifier()
    {
        refined_topo = 0;
        max_priority_queue_size = 0;
    }

    //function used to topologically simplify the model
    void exec_local_topological_simplification(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,
                                               OpType operation, int cache_size, coord_type persistence);
    coord_type get_average_persistence_value(IG &ig, Mesh &mesh);

    void exec_global_topological_simplification(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,
                                                OpType operation, int cache_size, coord_type persistence);
    
    inline void print_simplification_stats()
    {
        cerr<<"   performed simplification: "<<refined_topo<<endl;
        cerr<<"   maximum priority queue size: " <<max_priority_queue_size<< endl;
        refined_topo = max_priority_queue_size = 0;
    }

private:
    /// statistical variables for simplification
    utype refined_topo;
    utype max_priority_queue_size;

    /// ----- LOCAL TOPOLOGICAL SIMPLIFICATION FUNCTIONS ----- ///
    void local_topological_simplification(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division, Node_V &root, OpType operation, mig_cache &cache, coord_type persistence);
    void local_topological_simplification_leaf(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, IG& ig, local_VTstar_ET &local_rels, mig_cache &cache, Node_V &root, Spatial_Subdivision &division, coord_type persistence);
    void get_local_MIG(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, IG &ig, local_VTstar_ET &local_rels,
                       Node_V &root, Spatial_Subdivision &division, mig_cache &cache);
    void explore_desc1cell_local_mig(Node_V &n, const pair<itype,short> &v_pair, Mesh &mesh, Forman_Gradient &gradient, IG &ig, local_VTstar_ET &local_rels,
                                     iNode *saddle_node, Node_V &root, Spatial_Subdivision &division, mig_cache &cache);
    void explore_asc1cell_local_mig(Node_V &n, const pair<itype, short> &tf, iNode *saddle_node, itype first_t_id, leaf_ET &local_ef,
                                    Mesh &mesh, Forman_Gradient &gradient, IG &ig, Node_V &root, Spatial_Subdivision &division, mig_cache &cache);

    /// push into the priority queue only those arc below the avg
    void build_persistence_queue(priority_arcs_queue &q, IG &ig, Mesh &mesh, coord_type persistence);
    void simplify(priority_arcs_queue &queue, Node_V &n, Mesh &mesh, Forman_Gradient &gradient, IG& ig, local_VTstar_ET &local_rels, mig_cache &cache,
                  Node_V &root, Spatial_Subdivision &division, coord_type persistence);

    void contraction(nNode *extrema, iNode *saddle, priority_arcs_queue &q, IG &ig, Mesh &mesh, Forman_Gradient &gradient,
                     local_VTstar_ET &local_rels, mig_cache &cache, Node_V &n, Node_V &root, Spatial_Subdivision &division, coord_type persistence);
    void removal(nNode *extrema, iNode *saddle, priority_arcs_queue &q, IG &ig, Mesh &mesh, Forman_Gradient &gradient,
                 local_VTstar_ET &local_rels, mig_cache &cache, Node_V &n, Node_V &root, Spatial_Subdivision &division, coord_type persistence);

    void contraction_update_gradient(itype vertex, itype ex_minimum, itype next_vertex, ivect &critical_edge, Mesh &mesh, Forman_Gradient &gradient,
                                     local_VTstar_ET &local_rels, mig_cache &cache, Node_V &n, Node_V &root, Spatial_Subdivision &division);
    void removal_update_gradient(itype triangle, iNode *saddle, Mesh &mesh, Forman_Gradient &gradient, leaf_ET &local_ef, mig_cache &cache,
                                 Node_V &n, Node_V &root, Spatial_Subdivision &division);

    inline void remove_saddle_arcs(iNode *saddle, bool is_minimum, IG &ig)
    {
        set<Arc*> &arcs = saddle->getArcs(is_minimum);
        nNode* extreme_node;

        set<Arc*>::iterator it=arcs.begin();
        while(it!=arcs.end())
        {
            (*it)->setLabel(-1);
            if(is_minimum)
                extreme_node = ((nNode*)(*it)->getNode_i());
            else
                extreme_node = ((nNode*)(*it)->getNode_j());

            ig.removeArc(!is_minimum, (*it)); /// voglio zero se rimuovo un minimo, uno per il massimo
            extreme_node->removeArc(*it);
            saddle->removeArc(is_minimum, it); /// new -> elimino partendo dall'iteratore cosi' ho costo costante
            /// set again to begin
            it=arcs.begin();
        }
    }
    void remove_extreme_arcs(nNode* extrema, iNode* saddle, nNode *other_extrema, itype ending_path_simplex,
                             bool is_minimum, IG &ig, priority_arcs_queue &q, Mesh &mesh, coord_type persistence);
    /// return the position of the paired edge of the triangle - plus saves in the edge variable the paired edge
    inline int get_paired_edge(Forman_Gradient &gradient, Triangle& tri, itype t_id, ivect &edge)
    {
        for(int i=0; i<tri.vertices_num(); i++)
        {
            /// occorre controllare se l'edge e' accoppiato con il triangolo corrente
            tri.TE(i,edge);
            int t = gradient.convert_compressed_to_expand(t_id).get_edge_pair(tri.vertex_index(edge[0]),
                                                                     tri.vertex_index(edge[1]));
            if(t == 3)
            {
                return i;
            }
        }
        cout<<"[get_paired_edge] there is no paired edge into the triangle.. should not happen this"<<endl;
        cout<<"edge: "<<edge[0]<<" "<<edge[1]<<endl;
        cout<<"triangle: "<<t_id<<" "<<tri<<endl;
        int a; cin>>a;
        return -1;
    }


    /// ----- GLOBAL TOPOLOGICAL SIMPLIFICATION FUNCTIONS ----- ///
    void global_topological_simplification(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division, Node_V &root, OpType operation,
                                           mig_cache &cache, coord_type persistence);
    void global_topological_simplification_leaf(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, local_VTstar_ET &local_rels, mig_cache &cache,
                                                Node_V &root, Spatial_Subdivision &division, coord_type persistence);

    /// push into the priority queue only those arc below the avg
    void build_persistence_queue(priority_arcs_queue &q, leaf_ET &local_ef, Mesh &mesh, Forman_Gradient &gradient, coord_type persistence);
    void simplify(priority_arcs_queue &queue, Node_V &n, Mesh &mesh, Forman_Gradient &gradient, local_VTstar_ET &local_rels, mig_cache &cache,
                  Node_V &root, Spatial_Subdivision &division, coord_type persistence);
    void contraction(nNode *extrema, iNode *saddle, priority_arcs_queue &q, Mesh &mesh, Forman_Gradient &gradient,
                     local_VTstar_ET &local_rels, mig_cache &cache, Node_V &n, Node_V &root, Spatial_Subdivision &division, coord_type persistence);
    void removal(nNode *extrema, iNode *saddle, priority_arcs_queue &q, Mesh &mesh, Forman_Gradient &gradient,
                 local_VTstar_ET &local_rels, mig_cache &cache, Node_V &n, Node_V &root, Spatial_Subdivision &division, coord_type persistence);
    void remove_extreme_arcs(nNode* extrema, iNode* saddle, nNode *other_extrema, itype ending_path_simplex,
                             bool is_minimum, priority_arcs_queue &q, Node_V &n,Mesh &mesh, coord_type persistence);
};

#endif // FORMAN_GRADIENT_SIMPLIFICATION_H
