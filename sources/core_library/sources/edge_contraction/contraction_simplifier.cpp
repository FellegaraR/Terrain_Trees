/*
    This file is part of the Stellar library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The Stellar library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Stellar library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Stellar library.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "contraction_simplifier.h"
#include "terrain_trees/reindexer.h"
#include "utilities/container_utilities.h"



 void Contraction_Simplifier::find_candidate_edges(Node_V &n, Mesh &mesh,leaf_VT &vts, edge_queue& edges, contraction_parameters &params){
        map<ivect, coord_type> lengths;
        ivect e;
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {   
            RunIterator const& t_id = itPair.first;
            Triangle& t = mesh.get_triangle(*t_id);
            for(int i=0; i<t.vertices_num(); i++)
            {
                t.TE(i,e);
                if(n.indexes_vertex(e[1]))// e (v1,v2) is a candidate edge if at least v2 is in n
                {
                    map<ivect,coord_type>::iterator it = lengths.find(e);
                    if(it == lengths.end())
                    {
                    coord_type   length;
                    Vertex &v1=mesh.get_vertex(e[0]);
                    Vertex &v2=mesh.get_vertex(e[1]);
                    dvect dif = {v1.get_x()-v2.get_x(),v1.get_y()-v2.get_y(),v1.get_z()-v2.get_z()};
                    length = sqrt(dif[0]*dif[0]+dif[1]*dif[1]+dif[2]*dif[2]);
            //  Edge e((*it)[0],(*it)[1]);
                    lengths[e] = length;
                   //Edge edge_obj(e[0],e[1]);
                   if(length<params.get_maximum_length()){
                    Geom_Edge new_edge(e,length);
                    edges.push(&new_edge);
                    }
                    }
                }

            }   
        }
}


void Contraction_Simplifier::contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1,  Node_V &outer_v_block, edge_queue &edges,
                                           Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params)
{
    ivect et_vec;
    et_vec.push_back(et.first);
    et_vec.push_back(et.second);
    difference_of_vectors(vt0,et_vec); // vt0 now contains the difference VT0 - ET
    difference_of_vectors(vt1,et_vec); // vt1 now contains the difference VT1 - ET

    // contract v1 to v0.

    /// prior checking the d-1 faces we update
    /// (1) the corresponding vt0 relation (by adding the tops in vt1-et to vt0)
    /// (2) then the outer_v_block (if e is a cross edge)
    /// (3) and, finally, we add the new edges crossing b to the edge queue
    Contraction_Simplifier::update(e,vt0,vt1,n,outer_v_block,edges,mesh,params);

    // we remove v2 and the triangles in et
    Contraction_Simplifier::remove_from_mesh(e[1],et,mesh,params);  // What is e_top_id
    // finally we clear the VTop(v2)
    vt1.clear();
    //et.clear();
}



void Contraction_Simplifier::get_ET(ivect &e, ET &et, Node_V &n, Mesh &mesh, leaf_VT &vts)
{
    int other_v, local_v_id;
    if(n.indexes_vertex(e[0]))
    {
        other_v = e[1];
        local_v_id = e[0] - n.get_v_start();
    }
    else
    {
        other_v = e[0];
        local_v_id = e[1] - n.get_v_start();
    }

    VT &vt = vts[local_v_id];

    ivect et_tmp;
    
    for(unsigned i=0; i<vt.size();i++)
        {
            Triangle &t = mesh.get_triangle(vt[i]);
            if(t.has_vertex(other_v))
                et_tmp.push_back(vt[i]);
        }
    et=make_pair(et_tmp[0],et_tmp[1]);

}

void Contraction_Simplifier::clean_coboundary(VT &cob, Mesh &mesh)
{
  //  for(unsigned d=0; d<cob.size(); d++)
  //  {
        for(ivect_iter it=cob.begin(); it!=cob.end(); )
        {
            if(mesh.is_triangle_removed(*it))
                cob.erase(it);
            else
                ++it;
        }
  //  }
}

void Contraction_Simplifier::update(const ivect &e, VT& vt, VT& difference, Node_V &n, Node_V &v_block,
                                                            edge_queue &edges, Mesh &mesh, contraction_parameters &params)
{
    set<ivect> e_set; /// we insert the new edges first in this set to avoid duplicate insertions in the queue

    // for(unsigned d=0; d<difference.size(); d++)
    // {
        for(ivect_iter it=difference.begin(); it!=difference.end(); ++it)
        {
            Triangle &t = mesh.get_triangle(*it);

            /// before updating the top simplex, we check
            /// if the leaf block indexing e[0] does not contain the current top d-simplex we have to add it
            /// NOTA: there is one possible case.. as leaf block n already indexes the top simplices in e[1]
            /// NOTA2: we have just to check that n and v_block are different (as if they are equal the edge is internal in n)
            if(n.get_v_start() != v_block.get_v_start() && !v_block.indexes_triangle_vertices(t))
            {
//                if(debug)
//                    cout<<"[add top to leaf] "<<t<<" -> "<<*v_block<<endl;
                v_block.add_triangle(*it);
            }

            /// then we update the top simplex changing e[1] with e[0]
            int pos = t.vertex_index(e[1]);
            t.setTV(pos,e[0]);

            /// we have to add the new edges in the queue
            ivect new_e; new_e.assign(2,0);
            for(int i=0; i<t.vertices_num(); i++)
            {
                if(i!=pos)
                {
                    t.TE(i,new_e);  //t.TE(new_e,pos,i); need to check

                    if(n.indexes_vertex(new_e[1])) /// we process an edge only if it has all the extrema already processed "DON'T Understand"
                        e_set.insert(new_e);
                }
            }
        }
    //}

    /// we push the new "unique" edges in the queue
    for(auto it=e_set.begin(); it!=e_set.end(); ++it)
    {
        double length;
        Vertex &v1=mesh.get_vertex((*it)[0]);
        Vertex &v2=mesh.get_vertex((*it)[1]);
        dvect dif = {v1.get_x()-v2.get_x(),v1.get_y()-v2.get_y(),v1.get_z()-v2.get_z()};
        length = sqrt(dif[0]*dif[0]+dif[1]*dif[1]+dif[2]*dif[2]);
        ivect e{(*it)[0],(*it)[1]};
        Geom_Edge new_edge(e,length);
        edges.push(&new_edge);
    }

    /// finally we update the VT relation of e[0]
    unify_vectors(vt,difference);
}


void Contraction_Simplifier::remove_from_mesh(int to_delete_v,  ET &et, Mesh &mesh, contraction_parameters &params)
{
    mesh.remove_triangle(et.first);
    mesh.remove_triangle(et.second);

    mesh.remove_vertex(to_delete_v);
    params.increment_contracted_edges_counter();
}

 bool Contraction_Simplifier::link_condition(int v0, int v1, VT &vt0, VT &vt1,Mesh &mesh){

    iset vv0,vv1;
    for (int i=0;i<vt0.size();i++){
        Triangle t=mesh.get_triangle(vt0[i]);
    int v0_id=t.vertex_index(v0);
    vv0.insert(t.TV((v0_id+1)%3));
    vv0.insert(t.TV((v0_id+2)%3));
}
    for (int i=0;i<vt1.size();i++){
        Triangle t=mesh.get_triangle(vt1[i]);
        int v1_id=t.vertex_index(v1);
        vv1.insert(t.TV((v1_id+1)%3));
        vv1.insert(t.TV((v1_id+2)%3));
    }
    int counter=0;
    for(iset_iter it=vv1.begin();it!=vv1.end();it++){
        if(vv0.find(*it)!=vv0.end()){
            counter++;

        }

    }   

return counter<=2;
}