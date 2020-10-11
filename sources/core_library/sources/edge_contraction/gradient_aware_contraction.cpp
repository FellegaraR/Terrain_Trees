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


    cout<<"Simplification in leaf."<<endl;
boost::dynamic_bitset<> is_v_border(v_end-v_start);
//cout<<"Simplification in leaf."<<endl;
leaf_VT local_vts(v_range,VT());
n.get_VT_and_border(local_vts,is_v_border,mesh);
cout<<"Extracted VT and border edges"<<endl;
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

         //   cout<<"Vertex removed"<<endl;
         delete current;
        // if(edges_contracted_leaf>edge_num*0.2)
           // cout<<"Num of deleted in a leaf"<<edges_contracted_leaf<<"; 20% of the queue num:"<<edge_num*0.2<<endl;
            continue;

        }

        ivect sorted_e=e;
        std::sort(sorted_e.begin(),sorted_e.end());
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
        if(link_condition(e[0],e[1],*vt0,*vt1,et,mesh)&&valid_gradient_configuration(e[0],e[1],*vt0,*vt1,et,is_v_border[e[0]-v_start],is_v_border[e[1]-v_start],gradient,mesh)){
        contract_edge(e,et,*vt0,*vt1,*outer_v_block,edges,n,mesh,cache,params,gradient);
        edges_contracted_leaf++;
    // break;
        }
    // cout<<"Number of edges remaining:"<<edges.size()<<endl;
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

    cout<<"[debug]checking edge "<<v1<<", "<<v2<<endl;
    if(v1_is_border||v2_is_border){
      cout<<"border edge"<<endl;
    return false;}
    if(vt1.size()<4||vt2.size()<4){

        cout<<"less than 4 triangles"<<endl;
        return false;
    }
    
    int t1=et.first;
    int t2=et.second;
    if(gradient.is_triangle_critical(t1)||gradient.is_triangle_critical(t2))
    {
        cout<<"t1 or t2 is critical"<<endl;
        return false;
    }
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
        {
            cout<<"v3_sin:"<<v3_sin<<" v2:"<<v2<<endl;
            cout<<"v3 sin pair is v2"<<endl;  // why this case? 
      //      return false;
        }

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
        {
            cout<<"v3 des pair is v2"<<endl; // why this case? 
           // return false;
        }
    itype t3_des=-1;
    itype t3_sin=-1;
    bool v2_is_critical=true;
    bool edge0_is_critical=true;
    for(int i=0; i<vt2.size(); i++){
        
        Triangle &t = mesh.get_triangle(vt2[i]);
        int v2_pos=t.vertex_index(v2);
        for(int j=0;j<t.vertices_num();j++){
            if(j!=v2_pos)
                {
                    //vv2.insert(t.TV(j));
                t.TE(j,new_e);
                v2_edges.insert(new_e);            
                if(!gradient.is_edge_critical(new_e,vt2[i],mesh))
                edge0_is_critical = false;
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
            if(v2==71)
            cout<<"[DEBUG]"<<vt2[i]<<endl;
        if(gradient.is_triangle_critical(vt2[i])) {

             cout<<"vt2 is critical"<<endl;
            return false;}
        //Instead of searching for vtstar, we check all the triangles here
        if(!gradient.is_vertex_critical(v2,vt2[i],mesh)) 
        v2_is_critical=false;
    }

    if(v2_is_critical)
     {
         cout<<"v2 is critical"<<endl;
         return false;}
    if(edge0_is_critical)
    {
        cout<<"vv(v2) has critical edge"<<endl;
        return false;
    }

    ivect edge;
    edge={v1,v2};
    std::sort(edge.begin(),edge.end());
    if(gradient.is_edge_critical(edge,et,mesh))
    return false;
    itype t3_adj_sin=-1;
    itype t3_adj_des=-1;
    bool edge1_critical=true;
    bool edge2_critical=true;
    for(int i=0; i<vt1.size(); i++){
        if(gradient.is_triangle_critical(vt1[i])){
            cout<<"vt1 is critical"<<endl;
         return false;}
        for(int j=0; j<3; j++){
            int vid=mesh.get_triangle(vt1[i]).TV(j);
            if(vid==v3_des){
                if(vt1[i]!=t2)
                t3_adj_des=vt1[i];
                ivect edge1;
                edge1 ={v1,v3_des};
                sort(edge1.begin(),edge1.end());
            if(!gradient.is_edge_critical(edge1,vt1[i],mesh))
                edge1_critical=false;
            }
            else if(vid==v3_sin){
                if(vt1[i]!=t1)
                 t3_adj_sin=vt1[i];
                ivect edge1;
                edge1 ={v1,v3_sin};
                std::sort(edge1.begin(),edge1.end());
            if(!gradient.is_edge_critical(edge1,vt1[i],mesh))
                edge2_critical=false;
            }

        }
    }
    if(edge1_critical||edge2_critical)
     {
         cout<<"edge is critical"<<endl;
         return false;}
    //Check if v1 is point to v3_sin
    itype v1_pair=-1;
    itype v1_id=mesh.get_triangle(t1).vertex_index(v1);
    short v1_pair_id=gradient.convert_compressed_to_expand(t1).get_vertex_pair(v1_id);
    if(v1_pair_id==-1)
    {
        v1_id=mesh.get_triangle(t2).vertex_index(v1);
        v1_pair_id=gradient.convert_compressed_to_expand(t2).get_vertex_pair(v1_id);
        v1_pair= (v1_pair_id!=-1) ? mesh.get_triangle(t2).TV(v1_pair_id):-1;

    }
    else{
        v1_pair= mesh.get_triangle(t1).TV(v1_pair_id);
    }

    //Check if v2 is point to v3_sin/v3_des
    itype v2_pair=-1;
    itype v2_id=mesh.get_triangle(t1).vertex_index(v2);
    short v2_pair_id=gradient.convert_compressed_to_expand(t1).get_vertex_pair(v2_id);
    if(v2_pair_id==-1)
    {
        v2_id=mesh.get_triangle(t2).vertex_index(v2);
        v2_pair_id=gradient.convert_compressed_to_expand(t2).get_vertex_pair(v2_id);
        v2_pair= (v2_pair_id!=-1) ? mesh.get_triangle(t2).TV(v2_pair_id):-1;

    }
    else{
        v2_pair= mesh.get_triangle(t1).TV(v2_pair_id);
    }

    if(v1_pair!=v2&&v2_pair!=v1){
        cout<<"edge should be paired with either v1 or v2"<<endl;
        cout<<"v1: "<<v1<<" v1_pair:"<<v1_pair<<endl;
        cout<<"v2: "<<v2<<" v2_pair:"<<v2_pair<<endl;
        return false;
    }

    if(v3_sin_pair!=-1&& v3_sin_pair==v2){
        cout<<"v3_sin is paired with v2"<<endl;
        cout<<"t3_adj_sin:"<<t3_adj_sin<<endl;
     gradient.update_VE_adj_T(t3_adj_sin,v3_sin,v1,mesh,gradient);   
    }
    else if (v3_sin_pair!=-1&& v3_sin_pair==v1)
    {
           cout<<"v3_sin is paired with v1"<<endl;
           cout<<"v3_sin:"<<v3_sin<<" v3_sin_pair:"<<v3_sin_pair<<endl;
         cout<<"t3_sin:"<<t3_sin<<endl;
     gradient.update_VE_adj_T(t3_sin,v3_sin,v2,mesh,gradient);
    }
    else if (v1_pair!=-1&&v1_pair==v3_sin)
    {
        cout<<"v1 is paired with v3_sin"<<endl;
        cout<<"v1:"<<v1<<" v1_pair:"<<v1_pair<<endl;
     gradient.update_VE_adj_T(t3_sin,v2,v3_sin,mesh,gradient);
    }
    else if(v1_pair==v2&&v2_pair==v3_sin)
    {
        cout<<"v2 is paired with v3_sin"<<endl;
        cout<<"v2:"<<v2<<" v2_pair:"<<v2_pair<<endl;
        gradient.update_VE_adj_T(t3_adj_sin,v1,v3_sin,mesh,gradient);
    }

    if(v3_des_pair!=-1&& v3_des_pair==v2){
        // gradient.update_VE_adj_T(t3,v3_des,v1,mesh);
        cout<<"v3_des is paired with v2"<<endl;
         cout<<"t3_adj_des:"<<t3_adj_des<<endl;
     gradient.update_VE_adj_T(t3_adj_des,v3_des,v1,mesh,gradient);   
    }
    else if (v3_des_pair!=-1&& v3_des_pair==v1)
    {
        cout<<"v3_des is paired with v1"<<endl;
         cout<<"t3_des:"<<t3_des<<endl;
      //  gradient.update_VE_adj_T(t3,v3_des,v2,mesh);
     gradient.update_VE_adj_T(t3_des,v3_des,v2,mesh,gradient);
    }
    else if (v1_pair!=-1&&v1_pair==v3_des)
    {
         cout<<"v1 is paired with v3_des"<<endl;
        //gradient.update_VE_adj_T(t3,v2,v3_des,mesh);
        cout<<"v1:"<<v1<<" v1_pair:"<<v1_pair<<endl;
     gradient.update_VE_adj_T(t3_des,v2,v3_des,mesh,gradient);
    }
    else if(v1_pair==v2&&v2_pair==v3_des)
    {
        cout<<"v2 is paired with v3_des"<<endl;
        cout<<"v2:"<<v2<<" v2_pair:"<<v2_pair<<endl;
        gradient.update_VE_adj_T(t3_adj_des,v1,v3_des,mesh,gradient);
    }
    return true;
    //Triangle 


}



void Gradient_Aware_Simplifier::update(const ivect &e, VT& vt, VT& difference, Node_V &n, Node_V &v_block, edge_queue &edges,  Mesh &mesh, contraction_parameters &params)
{
   set<ivect> e_set; /// we insert the new edges first in this set to avoid duplicate insertions in the queue
    
    if(params.is_QEM()){
        initialQuadric[e[0]]+=initialQuadric[e[1]];
        set<ivect> v1_related_set;
        for(int i=0;i<vt.size();i++){
                ivect new_e; new_e.assign(2,0);
                Triangle &t = mesh.get_triangle(vt[i]);
                int v1_pos=t.vertex_index(e[0]);
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

            if((error-params.get_maximum_limit()<SMALL_TOLER)&&n.indexes_vertex(e[1])){
                    updated_edges[*it]=error;
              cout<<"["<<e[0]-1<<","<<e[1]-1<<"]  Error will be introduced: "<<error<<endl;

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

            if((error-params.get_maximum_limit()<SMALL_TOLER)&&n.indexes_vertex(e[1])){
                if(updated_edges.find(*it)!=updated_edges.end())
                {
                    updated_edges[*it]=error;
                }

            cout<<"["<<e[0]-1<<","<<e[1]-1<<"]  Error will be introduced: "<<error<<endl;

                 edges.push(new Geom_Edge(e,error));
        }
                 
       
        }
        else{
        dvect dif = {v1.get_x()-v2.get_x(),v1.get_y()-v2.get_y(),v1.get_z()-v2.get_z()};
        value = sqrt(dif[0]*dif[0]+dif[1]*dif[1]+dif[2]*dif[2]);
         e={(*it)[0],(*it)[1]};
         
        if((value-params.get_maximum_limit()<SMALL_TOLER)&&n.indexes_vertex(e[1])){
        // Geom_Edge new_edge(e,length);
 //    cout<<"["<<e[0]<<","<<e[1]<<"]  Edge length: "<<value<<endl;
        edges.push(new Geom_Edge(e,value));
        }
        
        }
    }

    /// finally we update the VT relation of e[0]
    unify_vectors(vt,difference);



}