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
            if(mesh.is_triangle_removed(*t_id)){
              //  cout<<"triangle removed"<<endl;
                continue;

            }
            Triangle& t = mesh.get_triangle(*t_id);
            
            for(int i=0; i<t.vertices_num(); i++)
            {
                t.TE(i,e);
                // if(mesh.get_vertex(e[0]).get_z()>mesh.get_vertex(e[1]).get_z()){
                //     int tmp=e[1];
                //     e[1]=e[0];
                //     e[0]=tmp;

                // }
              
                if(n.indexes_vertex(e[1]))// e (v1,v2) is a candidate edge if at least v2 is in n
                {
                    map<ivect,coord_type>::iterator it = lengths.find(e);
                   // cout<<e[0]<<" and "<<e[1]<<endl;
                    if(it == lengths.end())
                    {
                    
                    coord_type   length;
                    Vertex &v1=mesh.get_vertex(e[0]);
                    Vertex &v2=mesh.get_vertex(e[1]);
                    dvect dif = {v1.get_x()-v2.get_x(),v1.get_y()-v2.get_y(),v1.get_z()-v2.get_z()};
                  //  cout<<dif[0]<<", "<<dif[1]<<", "<<dif[2]<<endl;
                    length = sqrt(dif[0]*dif[0]+dif[1]*dif[1]+dif[2]*dif[2]);
                  cout<<"["<<e[0]<<","<<e[1]<<"]  Edge length: "<<length<<endl;
                    //  Edge e((*it)[0],(*it)[1]);
                    lengths[e] = length;
                   //Edge edge_obj(e[0],e[1]);
                   if(length<params.get_maximum_limit()){
                       
                    edges.push(new Geom_Edge(e,length));
                   //     cout<<"ENQUEUE"<<endl;
                    }
                    }
                }

            }   
        }
        cout<<"======NEXT NODE======"<<endl;
}


void Contraction_Simplifier::contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1,  Node_V &outer_v_block, edge_queue &edges,
                                           Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params)
{
    cout<<"[EDGE CONTRACTION] v1 and v2:"<<e[0]<<", "<<e[1]<<endl;
   // cout<<"[NOTICE] Contract Edge"<<endl;
    ivect et_vec;
    et_vec.push_back(et.first);
    if(et.second!=-1)
        et_vec.push_back(et.second);

    // cout<<"[DEBUG]ET:"<<endl;
    // if(et.first!=-1){
    //     Triangle t=mesh.get_triangle(et.first);

    //     cout<<"T1 "<<et.first<<": "<<t.TV(0)<<","<<t.TV(1)<<","<<t.TV(2)<<endl;
    //     }
    // if(et.second!=-1){
    //     Triangle t=mesh.get_triangle(et.second);

    //     cout<<"T2 "<<et.second<<": "<<t.TV(0)<<","<<t.TV(1)<<","<<t.TV(2)<<endl;
    //     }
    difference_of_vectors(vt0,et_vec); // vt0 now contains the difference VT0 - ET
    difference_of_vectors(vt1,et_vec); // vt1 now contains the difference VT1 - ET

   // cout<<"VT1 size: "<<vt0.size()<<" VT2 size: "<<vt1.size()<<endl;
    // contract v1 to v0.


    
    /// prior checking the d-1 faces we update
    /// (1) the corresponding vt0 relation (by adding the triangles in vt1-et to vt0)
    /// (2) then the outer_v_block (if e is a cross edge)
    /// (3) and, finally, we add the new edges crossing b to the edge queue
    Contraction_Simplifier::update(e,vt0,vt1,n,outer_v_block,edges,mesh,params);

    // we remove v2 and the triangles in et
    Contraction_Simplifier::remove_from_mesh(e[1],et,mesh,params); 
    // finally we clear the VT(v2)
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
   // cout<<"vt size:"<<vt.size()<<endl;

    for(unsigned i=0; i<vt.size();i++)
        {
            Triangle &t = mesh.get_triangle(vt[i]);
            if(t.has_vertex(other_v))
                et_tmp.push_back(vt[i]);
        }
       // cout<<"et size:"<<et_tmp.size()<<endl;
    if(et_tmp.size()==2)
    et=make_pair(et_tmp[0],et_tmp[1]);
    else if(et_tmp.size()==0)// Can be deleted, since the default is (-1,-1)
    {
        et=make_pair(-1,-1);
    }
    else
    {
        et=make_pair(et_tmp[0],-1);
    }
    
}

void Contraction_Simplifier::clean_coboundary(VT &cob, Mesh &mesh)
{

        for(ivect_iter it=cob.begin(); it!=cob.end(); )
        {
            if(mesh.is_triangle_removed(*it))
                cob.erase(it);
            else
                ++it;
        }

}

void Contraction_Simplifier::update(const ivect &e, VT& vt, VT& difference, Node_V &n, Node_V &v_block,
                                                            edge_queue &edges, Mesh &mesh, contraction_parameters &params)
{
    set<ivect> e_set; /// we insert the new edges first in this set to avoid duplicate insertions in the queue
    
    if(params.is_QEM()){
        initialQuadric[e[0]]+=initialQuadric[e[1]];

        set<ivect> v1_related_set;
        for(int i=0;i<vt.size();i++){

                ivect new_e; new_e.assign(2,0);
                Triangle &t = mesh.get_triangle(vt[i]);
                int v1_pos=t.vertex_index(e[0]);

          //      if(e[0]==24993||e[1]==24995)
                // {
                //     cout<<"[DEBUG]Triangle "<<vt[i]<<": "<<t.TV(0)<<","<<t.TV(1)<<", "<<t.TV(2)<<endl;
                //     }

                for(int i=0; i<t.vertices_num(); i++)
            {
                if(i!=v1_pos)
                {
                    t.TE(i,new_e);  //t.TE(new_e,pos,i); need to check

                    //cout<<"New edge"<<new_e[0]<<", "<<new_e[1]<<endl;
               //     if(n.indexes_vertex(new_e[1])) /// we process an edge only if it has all the extrema already processed 
                        v1_related_set.insert(new_e);
                       // updated_edges[new_e]=-1;
                }
            }
   
        }
        for(auto it=v1_related_set.begin(); it!=v1_related_set.end(); ++it)
    {
        //Calculate length
        double value;
        ivect e;
        Vertex &v1=mesh.get_vertex((*it)[0]);
        Vertex &v2=mesh.get_vertex((*it)[1]);
            int new_vertex_pos=-1;
              double error = compute_error((*it)[0],(*it)[1],mesh,new_vertex_pos);
          //    cout<<"[DEBUG] calculated error: "<<error<<endl;
              assert(new_vertex_pos!=-1);
                if(new_vertex_pos==1)
                {
                    e={(*it)[1],(*it)[0]};
                }
                else
                    e={(*it)[0],(*it)[1]};

            if(error<params.get_maximum_limit()&&n.indexes_vertex(e[1])){
                    updated_edges[*it]=error;
                cout<<"["<<e[0]<<","<<e[1]<<"]  Error will be introduced: "<<error<<endl;

                 edges.push(new Geom_Edge(e,error));
        }
                 
    }



    }

        for(ivect_iter it=difference.begin(); it!=difference.end(); ++it)
        {
            Triangle &t = mesh.get_triangle(*it);


            /// before updating the triangle, we check
            /// if the leaf block indexing e[0] does not contain the current triangle we have to add it
            /// NOTA: there is one possible case.. as leaf block n already indexes the triangle in e[1]
            /// NOTA2: we have just to check that n and v_block are different (as if they are equal the edge is internal in n)
            if(n.get_v_start() != v_block.get_v_start() && !v_block.indexes_triangle_vertices(t))
            {
//                if(debug)
//                    cout<<"[add top to leaf] "<<t<<" -> "<<*v_block<<endl;
                v_block.add_triangle(*it);
            }

            /// then we update the triangle changing e[1] with e[0]
            int pos = t.vertex_index(e[1]);
            t.setTV(pos,e[0]);
            dvect diff(4,0.0);
            if(params.is_QEM()){

            ivect new_e; new_e.assign(2,0);
            for(int i=0; i<t.vertices_num(); i++)
            {
                if(i!=pos)
                {
                    t.TE(i,new_e);  //t.TE(new_e,pos,i); need to check
                 //   cout<<"New edge"<<new_e[0]<<", "<<new_e[1]<<endl;
               //     if(n.indexes_vertex(new_e[1])) /// we process an edge only if it has all the extrema already processed 
                        e_set.insert(new_e);
                }
            }
            }
            else{
            /// we have to add the new edges in the queue
            ivect new_e; new_e.assign(2,0);
            for(int i=0; i<t.vertices_num(); i++)
            {
                if(i!=pos)
                {
                    t.TE(i,new_e);  //t.TE(new_e,pos,i); need to check
               //     if(n.indexes_vertex(new_e[1])) /// we process an edge only if it has all the extrema already processed 
                        e_set.insert(new_e);
                }
            }
            }
        }


    /// we push the new "unique" edges in the queue
    for(auto it=e_set.begin(); it!=e_set.end(); ++it)
    {

        //Calculate length
        double value;
        ivect e;
        Vertex &v1=mesh.get_vertex((*it)[0]);
        Vertex &v2=mesh.get_vertex((*it)[1]);
        if(params.is_QEM()){
            int new_vertex_pos=-1;
              double error = compute_error((*it)[0],(*it)[1],mesh,new_vertex_pos);
          //    cout<<"[DEBUG] calculated error: "<<error<<endl;
         
              assert(new_vertex_pos!=-1);
                if(new_vertex_pos==1)
                {
                    e={(*it)[1],(*it)[0]};
                }
                else
                    e={(*it)[0],(*it)[1]};

            if(error<params.get_maximum_limit()&&n.indexes_vertex(e[1])){
                if(updated_edges.find(*it)!=updated_edges.end())
                {
                    updated_edges[*it]=error;
                }

                cout<<"["<<e[0]<<","<<e[1]<<"]  Error will be introduced: "<<error<<endl;

                 edges.push(new Geom_Edge(e,error));
        }
                 
       
        }
        else{
        dvect dif = {v1.get_x()-v2.get_x(),v1.get_y()-v2.get_y(),v1.get_z()-v2.get_z()};
        value = sqrt(dif[0]*dif[0]+dif[1]*dif[1]+dif[2]*dif[2]);
         e={(*it)[0],(*it)[1]};
         
        if(value<params.get_maximum_limit()&&n.indexes_vertex(e[1])){
        // Geom_Edge new_edge(e,length);
        edges.push(new Geom_Edge(e,value));}
        
        }
    }

    /// finally we update the VT relation of e[0]
    unify_vectors(vt,difference);
}


void Contraction_Simplifier::remove_from_mesh(int to_delete_v,  ET &et, Mesh &mesh, contraction_parameters &params)
{
    if(et.first!=-1)
    {mesh.remove_triangle(et.first);
    params.increment_counter();
    }
   if(et.second!=-1){
    mesh.remove_triangle(et.second);
   params.increment_counter();
   }
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
  //  cout<<v1<<"'s VV size: "<<vv1.size()<<endl;
    for(iset_iter it=vv1.begin();it!=vv1.end();it++){
        if(vv0.find(*it)!=vv0.end()){
            counter++;

        }

    }   

return counter<=2;
}


void Contraction_Simplifier::update_new(const ivect &e, VT& vt, VT& difference, Node_V &n, Node_V &v_block, edge_queue &edges,
                                          Mesh &mesh, contraction_parameters &params, int new_vertex_pos){
    
            }

            
void Contraction_Simplifier::simplify(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli)
{   

    cerr<<"==Homology preserving simplification - weak-link condition=="<<endl;

    cerr<<"[NOTICED] Cache size: "<<cli.cache_size<<endl;
    LRU_Cache<int,leaf_VT> cache(cli.cache_size); // the key is v_start while the value are the VT relations
    contraction_parameters params;
    params.set_maximum_limit(cli.maximum_limit);
    if(cli.QEM_based)
        params.queue_criterion_QEM();
    else
    {
        params.queue_criterion_length();
    }
    
    Timer time;
    int simplification_round;
    int round = 1;
    if(params.is_QEM()){

         trianglePlane =vector<dvect>(mesh.get_triangles_num(),dvect(4,0));
        initialQuadric = vector<Matrix>(mesh.get_vertices_num()+1,Matrix(0.0));
        cout<<"=========Calculate triangle plane========"<<endl;
        compute_triangle_plane(mesh,trianglePlane);
        cout<<"=========Calculate initial QEM========"<<endl;
        compute_initial_QEM(mesh,trianglePlane);

        }
    time.start();
    while(1)
    {
        simplification_round = params.get_contracted_edges_num();  //checked edges

        /// HERE YOU NEED TO DEFINE A PROCEDURE FOR SIMPLIFY THE TIN BY USING THE SPATIAL INDEX
        this->simplify_compute(tree.get_root(),mesh,cache,tree.get_subdivision(),params,tree);

        cout<<"Num of edges enqueued:"<<params.get_sum_edge_queue_sizes()<<endl;
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
        //    if(round==2)
        //      break;
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

  void Contraction_Simplifier::simplify_compute(Node_V &n,  Mesh &mesh, LRU_Cache<int, leaf_VT> &cache,Spatial_Subdivision &division,  contraction_parameters &params, PRT_Tree &tree)
{
      if (n.is_leaf())
    {
        if(params.is_QEM()==true)
        simplify_leaf_QEM(n,mesh,cache,params,tree);
        else
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

void Contraction_Simplifier::simplify_leaf_QEM(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params,PRT_Tree& tree){

if(!n.indexes_vertices())
     return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();
    itype v_range = v_end - v_start;


    //cout<<"Simplification in leaf."<<endl;
    leaf_VT local_vts(v_range,VT());
    n.get_VT(local_vts,mesh);
    // Create a priority queue of candidate edges
    edge_queue edges;
    find_candidate_edges_QEM(n,mesh,local_vts,edges,params);
    int edge_num=edges.size();
    int edges_contracted_leaf=0;
    cout<<"Edge number:"<<edges.size()<<endl;
    params.add_edge_queue_size(edges.size());
    while(!edges.empty())
    {
        Geom_Edge* current = edges.top();
         ivect e=current->edge;
  //    cout<<"Start contraction."<<endl;
  //  cout<<"Edge Length:"<<current->val<<endl;

        edges.pop();

        if (mesh.is_vertex_removed(e[0])||mesh.is_vertex_removed(e[1])){


     //   cout<<"[DEBUG] edge not complete: "<<e[0]<<", "<<e[1]<<endl;
         //   cout<<"Vertex removed"<<endl;
         delete current;
        // if(edges_contracted_leaf>edge_num*0.2)
           // cout<<"Num of deleted in a leaf"<<edges_contracted_leaf<<"; 20% of the queue num:"<<edge_num*0.2<<endl;
            continue;

        }

        ivect sorted_e=e;
        sort(sorted_e.begin(),sorted_e.end());
        auto it= updated_edges.find(sorted_e);
        if(it!=updated_edges.end()){
        //int tmp=-1;
       // double error = compute_error(e[0],e[1],mesh,tmp);
        if(it->second!=current->val)
        {
         // cout<<"[DEBUG] edge: "<<sorted_e[0]<<", "<<sorted_e[1]<<"; updated error: "<<it->second<<"old error: "<<current->val<<endl;
            delete current;
            continue;
        }
        }
        

        ET et(-1,-1);
        VT *vt0=NULL,*vt1=NULL;
        Node_V *outer_v_block=NULL;

        get_edge_relations(e,et,vt0,vt1,outer_v_block,n,mesh,local_vts,cache,params,tree);
        if(link_condition(e[0],e[1],*vt0,*vt1,mesh)){
        contract_edge(e,et,*vt0,*vt1,*outer_v_block,edges,n,mesh,cache,params);
        edges_contracted_leaf++;
    // break;
        }

    }



}


void Contraction_Simplifier::simplify_leaf(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params,PRT_Tree& tree){

if(!n.indexes_vertices())
     return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();
    itype v_range = v_end - v_start;


//cout<<"Simplification in leaf."<<endl;
leaf_VT local_vts(v_range,VT());
n.get_VT(local_vts,mesh);
// Create a priority quue of candidate edges
edge_queue edges;
find_candidate_edges(n,mesh,local_vts,edges,params);
int edge_num=edges.size();
int edges_contracted_leaf=0;
cout<<"Edge number:"<<edges.size()<<endl;
params.add_edge_queue_size(edges.size());
    while(!edges.empty())
    {
        Geom_Edge* current = edges.top();
         ivect e=current->edge;
  //    cout<<"Start contraction."<<endl;
  //  cout<<"Edge Length:"<<current->val<<endl;

        edges.pop();

        if (mesh.is_vertex_removed(e[0])||mesh.is_vertex_removed(e[1])){

         //   cout<<"Vertex removed"<<endl;
         delete current;
        // if(edges_contracted_leaf>edge_num*0.2)
           // cout<<"Num of deleted in a leaf"<<edges_contracted_leaf<<"; 20% of the queue num:"<<edge_num*0.2<<endl;
            continue;

        }

        ET et(-1,-1);
        VT *vt0=NULL,*vt1=NULL;
        Node_V *outer_v_block=NULL;

        get_edge_relations(e,et,vt0,vt1,outer_v_block,n,mesh,local_vts,cache,params,tree);
        if(link_condition(e[0],e[1],*vt0,*vt1,mesh)){
        contract_edge(e,et,*vt0,*vt1,*outer_v_block,edges,n,mesh,cache,params);
        edges_contracted_leaf++;
    // break;
        }

    }


// leaf_VV vvs;
// n.get_VV(vvs,mesh);
// for()


}

void Contraction_Simplifier::get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1, Node_V *&outer_v_block,
                                                Node_V &n, Mesh &mesh, leaf_VT &vts, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params,PRT_Tree &tree)
{

    //cout<<"[NOTICE]get edge relation"<<endl;
    outer_v_block = NULL;
    /// inverted order as I only need the block indexing v1
    vt1 = Contraction_Simplifier::get_VT(e[1],n,mesh,vts,cache,tree,outer_v_block,params);
    vt0 = Contraction_Simplifier::get_VT(e[0],n,mesh,vts,cache,tree,outer_v_block,params);



    Contraction_Simplifier::get_ET(e,et,n,mesh,vts);


}

 VT* Contraction_Simplifier::get_VT(int v_id, Node_V &n, Mesh &mesh, leaf_VT &vts, LRU_Cache<int,leaf_VT> &cache,
                        PRT_Tree &tree, Node_V *& v_block, contraction_parameters &params)
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
        }
            Contraction_Simplifier::clean_coboundary((it_c->second)[local_index],mesh);
            if(debug)
            cout<<"[NOTICE]cleaned coboundary"<<endl;
            // if(debug/* || v_id ==2355*/)
            //     cout<<"num_elem_in_vt: "<<get_num_elements_in_container_of_containers((it_c->second)[local_index])<<endl;
//                print_container_of_containers_content("CleanedVTop(2355) ",(it_c->second)[local_index]);
        
        return &(it_c->second)[local_index];
    }
}


void Contraction_Simplifier::update_mesh_and_tree(PRT_Tree &tree, Mesh &mesh, contraction_parameters &params)
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
    cout<<"number of surviving vertices:"<<surviving_vertices.size()<<endl;
    mu.clean_vertices_array(mesh,new_v_positions,surviving_vertices);
   cout<<"number of deleted triangles:"<< params.get_counter()<<endl;
    /// NEW: the update_and_compact procedure check internally if we have removed all the top d-simplices
    bool all_deleted = mu.update_and_clean_triangles_arrays(mesh,new_v_positions,new_t_positions,params.get_counter());
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


void Contraction_Simplifier::find_candidate_edges_QEM(Node_V &n, Mesh &mesh, leaf_VT &vts,edge_queue& edges, contraction_parameters &params){

      map<ivect, coord_type> edge_map;
        ivect e;
        int t_count=0;
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {   
            RunIterator const& t_id = itPair.first;
            if(mesh.is_triangle_removed(*t_id)){
              //  cout<<"triangle removed"<<endl;
                continue;

            }
            t_count++;
            Triangle& t = mesh.get_triangle(*t_id);
            
            for(int i=0; i<t.vertices_num(); i++)
            {
                t.TE(i,e);
 

              int new_vertex_pos=-1;
              double error = compute_error(e[0],e[1],mesh,new_vertex_pos);
              assert(new_vertex_pos!=-1);
                if(new_vertex_pos==1)
                {
                    int tmp=e[1];
                    e[1]=e[0];
                    e[0]=tmp;
                }
       
                if(n.indexes_vertex(e[1]))// e (v1,v2) is a candidate edge if at least v2 is in n
                {
                    map<ivect,coord_type>::iterator it = edge_map.find(e);
                   // cout<<e[0]<<" and "<<e[1]<<endl;
                    if(it == edge_map.end())
                    {
                    
                                   //  Edge e((*it)[0],(*it)[1]);
                    edge_map[e] = error;
                   //Edge edge_obj(e[0],e[1]);
                   if(error<params.get_maximum_limit()){
                    cout<<"["<<e[0]<<","<<e[1]<<"]  Error will be introduced: "<<error<<endl;

                    edges.push(new Geom_Edge(e,error));
                   //     cout<<"ENQUEUE"<<endl;
                    }
                    }
                }

            }   
        }
        cout<<"**** [Number] "<<edges.size()<<" edges enqueued. Start simplification.****"<<endl;
         cout<<"number of remaining triangles:"<<t_count<<endl;

        cout<<"======NEXT NODE======"<<endl;
}

void Contraction_Simplifier::compute_initial_QEM( Mesh &mesh, vector<dvect >& planes){

 for (int i = 1; i <= mesh.get_triangles_num(); i++)
    {
        /* faces are triangles */
        for (int j = 0; j < 3; j++)
        {
            double* a = &(planes[i-1][0]);
            initialQuadric[ mesh.get_triangle(i).TV(j) ] += Matrix(a);
        }
    }
}



void Contraction_Simplifier::compute_triangle_plane(Mesh &mesh,vector<dvect>& trPl){
double coords[3][3];
for(int i=1; i<=mesh.get_triangles_num(); i++){
        
        for(int v=0; v<3; v++){
            Vertex v1= mesh.get_vertex(mesh.get_triangle(i).TV(v));
            coords[0][v] = v1.get_x();
            coords[1][v] = v1.get_y();
            coords[2][v] = v1.get_z();

        }

        double a,b,c,m;

        a = (coords[1][1] - coords[1][0]) * (coords[2][2] - coords[2][0]) - (coords[2][1] - coords[2][0]) * (coords[1][2] - coords[1][0]);

        b = (coords[2][1] - coords[2][0]) * (coords[0][2] - coords[0][0]) - (coords[0][1] - coords[0][0]) * (coords[2][2] - coords[2][0]);

        c = (coords[0][1] - coords[0][0]) * (coords[1][2] - coords[1][0]) - (coords[1][1] - coords[1][0]) * (coords[0][2] - coords[0][0]);

        m = sqrt(a*a + b*b + c*c);
        a = a/m;
        b = b/m;
        c = c/m;

        trPl[i-1][0]=a;
        trPl[i-1][1]=b;
        trPl[i-1][2]=c;
        trPl[i-1][3]= -1*(a*coords[0][0] + b*coords[1][0] + c*coords[2][0]);
    }
}



double Contraction_Simplifier::compute_error(int v1, int v2, Mesh &mesh,int& new_vertex_pos){
 double min_error;
    Matrix q_bar;
    Matrix q_delta;


    /* computer quadric of virtual vertex vf */
    //cout<<"v1: "<<v1<<" and v2: "<<v2<<endl;
    q_bar = initialQuadric[v1] + initialQuadric[v2];


    q_delta = Matrix( q_bar[0], q_bar[1],  q_bar[2],  q_bar[3],
                      q_bar[4], q_bar[5],  q_bar[6],  q_bar[7],
                      q_bar[8], q_bar[9], q_bar[10], q_bar[11],
                             0,        0,	      0,        1);

    Vertex vertex_1 = mesh.get_vertex(v1);
    double vx1 = vertex_1.get_x();
    double vy1 = vertex_1.get_y();
    double vz1 = vertex_1.get_z();

    Vertex vertex_2 = mesh.get_vertex(v2);
    double vx2 = vertex_2.get_x();
    double vy2 = vertex_2.get_y();
    double vz2 = vertex_2.get_z();


        double vx3 = double (vx1+vx2)/2.0;
        double vy3 = double (vy1+vy2)/2.0;
        double vz3 = double (vz1+vz2)/2.0;

        double error1 = vertex_error(q_bar, vx1, vy1, vz1);
        double error2 = vertex_error(q_bar, vx2, vy2, vz2);
 

        min_error = std::min(error1,error2);
      
        if (error1 == min_error) {new_vertex_pos=0; min_error=error1; }
        else {new_vertex_pos=1; min_error=error2; }
 


  //  min_error = vertex_error(q_bar, new_vertex[0], new_vertex[1], new_vertex[2]);
if(min_error<=0.00000001)
    min_error=0;
    return min_error;

}