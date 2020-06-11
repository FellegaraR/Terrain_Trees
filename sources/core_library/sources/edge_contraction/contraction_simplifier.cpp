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
        this->simplify_compute(tree.get_root(),mesh,cache,tree.get_subdivision(),params);
        // PARTIAL SIMPLIFICATION STATS
        cerr<<"=== end-of-round "<<round<<") --> contracted edges: ";
        cerr<<params.get_contracted_edges_num()-simplification_round<<endl;
        round++;

        if(cli.debug_mode)
        {
            time.stop();
            time.print_elapsed_time("   [TIME] executing a simplification round: ");
            time.start();
            cerr << "   [RAM] peak for executing a simplification round: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
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

    cerr << "[RAM peak] for contracting a simplicial complex: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

    /// finally we have to update/compress the mesh and the tree
    Contraction_Simplifier::update_mesh_and_tree(tree,mesh,params);
}

 void simplify_compute(Node_V &n,  Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, Spatial_Subdivision &division, contraction_parameters &params)
{
      if (n.is_leaf())
    {
        simplify_leaf(n,mesh,cache,params);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                simplify_compute(*n.get_son(i),mesh,cache,division,params);
            }
        }
    }

}



static void simplify_leaf(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params){

if(!n.indexes_vertices())
     return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();
    itype v_range = v_end - v_start;

leaf_VT local_vts(v_range,VT());
n.get_VT(local_vts,mesh);
// Create a priority quue of candidate edges

    
}


static void find_candidate_edges(Node_V &n, Mesh &mesh, edge_queue& edges, contraction_parameters &params){

    
}


void Contraction_Simplifier::contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1, leaf_VT &vt, Node_V &outer_v_block, edge_queue &edges,
                                           Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params)
{
    ivect et_vec;
    et_vec.push_back(et.first);
    et_vec.push_back(et.second);
    difference_of_container_of_containers(vt0,et_vec); // vt0 now contains the difference VT0 - ET
    difference_of_container_of_containers(vt1,et_vec); // vt1 now contains the difference VT1 - ET

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

template<class T> void Contraction_Simplifier::get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1, Node_V *&outer_v_block,
                                                Node_V &n, Mesh &mesh, leaf_VT &vt, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params,T &tree)
{


    outer_v_block = NULL;
    /// inverted order as I only need the block indexing v1
    vt1 = Contraction_Simplifier::get_VT(e[1],n,mesh,vt,cache,tree,outer_v_block,params);
    vt0 = Contraction_Simplifier::get_VT(e[0],n,mesh,vt,cache,tree,outer_v_block,params);



    Contraction_Simplifier::get_ET(e,et,n,mesh,vt);


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
            cout<<"[get_VTop] "<<v_id<<" -> EXTERNAL VERTEX "<<n<<endl;

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

            if(debug/* || v_id ==2355*/)
                cout<<"num_elem_in_vt: "<<get_num_elements_in_container_of_containers((it_c->second)[local_index])<<endl;
//                print_container_of_containers_content("VTop(2355) ",(it_c->second)[local_index]);

            Contraction_Simplifier::clean_coboundary((it_c->second)[local_index],mesh);

            if(debug/* || v_id ==2355*/)
                cout<<"num_elem_in_vt: "<<get_num_elements_in_container_of_containers((it_c->second)[local_index])<<endl;
//                print_container_of_containers_content("CleanedVTop(2355) ",(it_c->second)[local_index]);
        }

        return &(it_c->second)[local_index];
    }
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
        Edge e((*it)[0],(*it)[1]);
        Geom_Edge new_edge(&e,length);
        edges.push(&new_edge);
    }

    /// finally we update the VT relation of e[0]
    unify_container_of_containers(vt,difference);
}


void Contraction_Simplifier::remove_from_mesh(int to_delete_v,  ET &et, Mesh &mesh, contraction_parameters &params)
{
    mesh.remove_triangle(et.first);
    mesh.remove_triangle(et.second);

    mesh.remove_vertex(to_delete_v);
    params.increment_contracted_edges_counter();
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
    boost::dynamic_bitset<> all_deleted = mu.update_and_clean_triangles_arrays(mesh,new_v_positions,new_top_positions,params.get_counters());
    time.stop();
    time.print_elapsed_time("[TIME] Compact and update mesh: ");

    cerr<<"[STAT] mesh "<<endl;
    cerr<<"------------"<<endl;
    mesh.print_mesh_stats(cerr);
    cerr<<"------------"<<endl;

    time.start();
//    cerr<<"[TREE] update indices in the tree"<<endl;
    tree.update_tree(tree.get_root(),new_v_positions,new_top_positions,all_deleted);
    time.stop();
    time.print_elapsed_time("[TIME] Update tree (top-simplices): ");

//    Reindexer r;
//    r.reorganize_index_and_mesh(tree,mesh,false);

    cerr << "[RAM peak] for updating the mesh and the tree: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
}
