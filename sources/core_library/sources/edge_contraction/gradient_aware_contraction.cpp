#include "gradient_aware_contraction.h"

void Gradient_Aware_Simplifier::gradient_aware_simplify(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli,Forman_Gradient &gradient)   
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
        this->simplify_compute(tree.get_root(),mesh,cache,tree.get_subdivision(),params,tree,gradient);

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

    }
    time.stop();
    if(!cli.debug_mode)
        time.print_elapsed_time("[TIME] Edge contraction simplification: ");

  //  cerr << "[RAM peak] for contracting a simplicial complex: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

    /// finally we have to update/compress the mesh and the tree
    Contraction_Simplifier::update_mesh_and_tree(tree,mesh,params);
}

  void Gradient_Aware_Simplifier::simplify_compute(Node_V &n,  Mesh &mesh, LRU_Cache<int, leaf_VT> &cache,Spatial_Subdivision &division,  contraction_parameters &params, PRT_Tree &tree,Forman_Gradient &gradient)
{
      if (n.is_leaf())
    {
        if(params.is_QEM()==true)
        simplify_leaf_QEM(n,mesh,cache,params,tree,gradient);
        else
        simplify_leaf(n,mesh,cache,params,tree,gradient);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                simplify_compute(*n.get_son(i),mesh,cache,division,params,tree,gradient);
            }
        }
    }

}

void Gradient_Aware_Simplifier::simplify_leaf_QEM(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params,PRT_Tree& tree,Forman_Gradient &gradient){

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
  //  cout<<"Edge error:"<<current->val<<endl;

        edges.pop();

        if (mesh.is_vertex_removed(e[0])||mesh.is_vertex_removed(e[1])){

     //       cout<<"skip current edge."<<endl;
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
        if(fabs(it->second-current->val)>SMALL_TOLER)
        {
           // cout<<"skip current edge."<<endl;
         // cout<<"[DEBUG] edge: "<<sorted_e[0]<<", "<<sorted_e[1]<<"; updated error: "<<it->second<<"old error: "<<current->val<<endl;
            delete current;
            continue;
        }
        }
        

        ET et(-1,-1);
        VT *vt0=NULL,*vt1=NULL;
        Node_V *outer_v_block=NULL;

        get_edge_relations(e,et,vt0,vt1,outer_v_block,n,mesh,local_vts,cache,params,tree);
        if(link_condition(e[0],e[1],*vt0,*vt1,et,mesh)){
        contract_edge(e,et,*vt0,*vt1,*outer_v_block,edges,n,mesh,cache,params,gradient);
        edges_contracted_leaf++;
    // break;
        }
     cout<<"Number of edges remaining:"<<edges.size()<<endl;
    }



}


void Gradient_Aware_Simplifier::simplify_leaf(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params,PRT_Tree& tree,Forman_Gradient &gradient){

if(!n.indexes_vertices())
     return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();
    itype v_range = v_end - v_start;

boost::dynamic_bitset<> is_v_border(v_end-v_start);
//cout<<"Simplification in leaf."<<endl;
leaf_VT local_vts(v_range,VT());
n.get_VT_and_border(local_vts,is_v_border,mesh);

// local_VTstar_ET all_rels;
// Forman_Gradient_Topological_Relations::get_VTstar_ET(all_rels,n,mesh,gradient);

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
        // cout<<"skip current edge"<<endl;
         delete current;
        // if(edges_contracted_leaf>edge_num*0.2)
           // cout<<"Num of deleted in a leaf"<<edges_contracted_leaf<<"; 20% of the queue num:"<<edge_num*0.2<<endl;
            continue;

        }

        ET et(-1,-1);
        VT *vt0=NULL,*vt1=NULL;
        Node_V *outer_v_block=NULL;
        
        get_edge_relations(e,et,vt0,vt1,outer_v_block,n,mesh,local_vts,cache,params,tree);
        if(link_condition(e[0],e[1],*vt0,*vt1,et,mesh)){
        contract_edge(e,et,*vt0,*vt1,*outer_v_block,edges,n,mesh,cache,params,gradient);
        edges_contracted_leaf++;
    // break;
        }

    }

}

void Gradient_Aware_Simplifier::contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1,  Node_V &outer_v_block, edge_queue &edges,
                                           Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params,Forman_Gradient &gradient)
{
    cout<<"[EDGE CONTRACTION] v1 and v2:"<<e[0]-1<<", "<<e[1]-1<<endl;
   // cout<<"[NOTICE] Contract Edge"<<endl;
    ivect et_vec;
    et_vec.push_back(et.first);
    if(et.second!=-1)
        et_vec.push_back(et.second);

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


bool Gradient_Aware_Simplifier::valid_gradient_configuration(int v1,int v2, VT &vt1, VT &vt2,ET& et ,bool v1_is_border, bool v2_is_border, Forman_Gradient &gradient, Mesh &mesh){

    if(v1_is_border||v2_is_border)
    return false;
    if(vt1.size()<4||vt2.size()<4)
    return false;
    int t1=et.first;
    int t2=et.second;
    if(gradient.is_triangle_critical(t1)||gradient.is_triangle_critical(t2))
    return false;

    int v3_sin, v3_des;
    v3_sin = v3_des = -1;
 //   iset vv2;
      ivect new_e; new_e.assign(2,0);
    set<ivect> v2_edges;
    



    short v3_sin_pair_id;
    for(int i=0; i<3; i++){
        if(mesh.get_triangle(t1).TV(i) != v1 && mesh.get_triangle(t1).TV(i) != v2){
            v3_sin = mesh.get_triangle(t1).TV(i);
            v3_sin_pair_id=gradient.convert_compressed_to_expand(t1).get_vertex_pair(i);
            break;
        }
    }

    itype v3_sin_pair= (v3_sin_pair_id!=-1) ? mesh.get_triangle(t1).TV(v3_sin_pair_id):-1;
    if(v3_sin_pair==v2)
        return false;

    short v3_des_pair_id;
    for(int i=0; i<3; i++){
        if(mesh.get_triangle(t2).TV(i) != v1 && mesh.get_triangle(t2).TV(i) != v2){
            v3_des = mesh.get_triangle(t2).TV(i);// i is vertex_index of v3_des in t2
            // we can then use it for checking gradient.
            v3_des_pair_id=gradient.convert_compressed_to_expand(t2).get_vertex_pair(i);
            break;
        }
    }
    itype v3_des_pair= (v3_des_pair_id!=-1) ? mesh.get_triangle(t2).TV(v3_des_pair_id):-1;
    if(v3_des_pair==v2)
        return false;

    itype t3_des=-1;
    itype t3_sin=-1;
    for(int i=0; i<vt2.size(); i++){
        
        Triangle &t = mesh.get_triangle(vt2[i]);
        int v2_pos=t.vertex_index(v2);
        for(int j=0;j<t.vertices_num();j++){
            if(j!=v2_pos)
                {
                    //vv2.insert(t.TV(j));
                t.TE(j,new_e);
                v2_edges.insert(new_e);            
                if(gradient.is_edge_critical(new_e,vt2[i],mesh))
                return false;
            if(t.TV(j)==v3_des&&vt2[i]!=t2)
            {
                t3_des=vt2[i];
            }
            else if(t.TV(j)==v3_sin&&vt2[i]!=t1)
            {
                t3_sin=vt2[i];
            }
                }
        }

        if(gradient.is_triangle_critical(vt2[i])) return false;
        //Instead of searching for vtstar, we check all the triangles here
        if(gradient.is_vertex_critical(v2,vt2[i],mesh));
    }



    ivect edge;
    edge={v1,v2};
    if(gradient.is_edge_critical(edge,et,mesh))
    return false;
    itype t3_adj_sin=-1;
    itype t3_adj_des=-1;

    for(int i=0; i<vt1.size(); i++){
        if(gradient.is_triangle_critical(vt1[i]))
         return false;
        for(int j=0; j<3; j++){
            int vid=mesh.get_triangle(vt1[i]).TV(j);
            if(vid==v3_des){
                if(vt1[i]!=t2)
                t3_adj_des=vt1[i];
                ivect edge;
                edge ={v1,v3_des};
            if(gradient.is_edge_critical(edge,vt1[i],mesh))
                return false;
            }
            else if(vid==v3_sin){
                if(vt1[i]!=t1)
                 t3_adj_sin=vt1[i];
                ivect edge;
                edge ={v1,v3_sin};
            if(gradient.is_edge_critical(edge,vt1[i],mesh))
                return false;
            }

        }
    }

    //Check if v1 is point to v3_sin
    itype v1_id=mesh.get_triangle(t1).vertex_index(v1);
    short v1_pair_id=gradient.convert_compressed_to_expand(t1).get_vertex_pair(v1_id);
    itype v1_pair= (v1_pair_id!=-1) ? mesh.get_triangle(t1).TV(v1_pair_id):-1;
    


    if(v3_sin_pair!=-1&& v3_sin_pair==v2){
     gradient.update_VE_adj_T(t3_adj_sin,v3_sin,v1,mesh,gradient);   
    }
    else if (v3_sin_pair!=-1&& v3_sin_pair==v1)
    {
     gradient.update_VE_adj_T(t3_sin,v3_sin,v2,mesh,gradient);
    }
    else if (v1_pair!=-1&&v1_pair==v3_sin)
    {
     gradient.update_VE_adj_T(t3_sin,v2,v3_sin,mesh,gradient);
    }


    if(v3_des_pair!=-1&& v3_des_pair==v2){
        // gradient.update_VE_adj_T(t3,v3_des,v1,mesh);
     gradient.update_VE_adj_T(t3_adj_des,v3_des,v1,mesh,gradient);   
    }
    else if (v3_des_pair!=-1&& v3_des_pair==v1)
    {
      //  gradient.update_VE_adj_T(t3,v3_des,v2,mesh);
     gradient.update_VE_adj_T(t3_des,v3_des,v2,mesh,gradient);
    }
    else if (v1_pair!=-1&&v1_pair==v3_des)
    {
        //gradient.update_VE_adj_T(t3,v2,v3_des,mesh);
     gradient.update_VE_adj_T(t3_des,v2,v3_des,mesh,gradient);
    }
    //Triangle 

    /* ******IN IA***********
    if(getVE(v3_sin) != NULL && getVE(v3_sin)->EV(1) == v2) return false;
    if(getVE(v3_des) != NULL && getVE(v3_des)->EV(1) == v2) return false;

    if(is_vertex_critical(v2))  //need vtstar of v2
        return false;

    Edge* edge1=getVE(v3_sin);
    Edge* edge2=getVE(v1);

    
    edge1=getVE(v3_des);
    	if(edge1 != NULL && *edge1==Edge(v3_des,v2)){

        int t3 = mesh->getTopSimplex(t2).TT(mesh->getTopSimplex(t2).vertex_index(v2));
        gradient.update_VE_adj_T(t3,v3_des,v1,mesh);
     
    }
    else if(edge1 != NULL && *edge1==Edge(v3_des,v1)){
        int t3 = mesh->getTopSimplex(t2).TT(mesh->getTopSimplex(t2).vertex_index(v1));
		gradient.update_VE_adj_T(t3,v3_des,v2,mesh);
               

    }
    else if(edge2 != NULL && *edge2==Edge(v3_des,v1)){
        int t3 = mesh->getTopSimplex(t2).TT(mesh->getTopSimplex(t2).vertex_index(v1));
		gradient.update_VE_adj_T(t3,v2,v3_des,mesh);

    }

    delete edge1;
    delete edge2;


    */

}



