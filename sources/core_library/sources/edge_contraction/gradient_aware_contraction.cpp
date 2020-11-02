#include "gradient_aware_contraction.h"

void Gradient_Aware_Simplifier::gradient_aware_simplify(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli,Forman_Gradient &gradient)   
{

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


 //   cout<<"Simplification in leaf."<<endl;
boost::dynamic_bitset<> is_v_border(v_end-v_start);
//cout<<"Simplification in leaf."<<endl;
leaf_VT local_vts(v_range,VT());
n.get_VT_and_border(local_vts,is_v_border,mesh);
//cout<<"Extracted VT and border edges"<<endl;
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
        if(fabs(it->second-current->val)>Zero)
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
        bool v1_is_border=false, v2_is_border=false;

        get_edge_relations(e,et,vt0,vt1,v1_is_border,v2_is_border,outer_v_block,n,mesh,local_vts,is_v_border,cache,params,tree);
        if(link_condition(e[0],e[1],*vt0,*vt1,et,mesh)&&valid_gradient_configuration(e[0],e[1],*vt0,*vt1,et,v1_is_border,v2_is_border,gradient,mesh)){
        contract_edge(e,et,*vt0,*vt1,*outer_v_block,edges,n,mesh,cache,params,gradient);
        edges_contracted_leaf++;
    // break;
        }
    // cout<<"Number of edges remaining:"<<edges.size()<<endl;
    }

// if(cache.find(v_start) != cache.end()){
//     cache.update(v_start,local_vts);
// }

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
//cout<<"Edge number:"<<edges.size()<<endl;
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
        bool v1_is_border=false, v2_is_border=false;
        
        get_edge_relations(e,et,vt0,vt1,v1_is_border,v2_is_border,outer_v_block,n,mesh,local_vts,is_v_border,cache,params,tree);
        if(link_condition(e[0],e[1],*vt0,*vt1,et,mesh)&&valid_gradient_configuration(e[0],e[1],*vt0,*vt1,et,v1_is_border,v2_is_border,gradient,mesh)){
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

    bool debug=false;
    cout<<"[debug]checking edge "<<v1<<", "<<v2<<endl;
    if(v1_is_border||v2_is_border){
      //cout<<"border edge"<<endl;
    return false;}
    if(vt1.size()<4||vt2.size()<4){
       // cout<<"less than 4 triangles"<<endl;
        return false;
    }
    
    int t1=et.first;
    int t2=et.second;
    if(gradient.is_triangle_critical(t1)||gradient.is_triangle_critical(t2))
    {
     //   cout<<"t1 or t2 is critical"<<endl;
        return false;
    }
    int v3_sin, v3_des;
    v3_sin = v3_des = -1;
 //   iset vv2;
    ivect new_e; new_e.assign(2,0);
   //set<ivect> v2_edges;
    
    if(debug){
        cout<<"[DEBUG]vt1 v1:"<<v1<<endl;
        for(auto it=vt1.begin();it!=vt1.end();it++){
            Triangle t=mesh.get_triangle(*it);
            cout<<t<<endl;

        }
        cout<<endl;

        cout<<"[DEBUG]vt2 v2:"<<v2<<endl;
        for(auto it=vt2.begin();it!=vt2.end();it++){
             Triangle t=mesh.get_triangle(*it);
            cout<<t<<endl;

        }
        cout<<endl;
    }


    short v3_sin_pair_id;
    for(int i=0; i<3; i++){
        if(mesh.get_triangle(t1).TV(i) != v1 && mesh.get_triangle(t1).TV(i) != v2){
            v3_sin = mesh.get_triangle(t1).TV(i);
            v3_sin_pair_id=gradient.convert_compressed_to_expand(t1).get_vertex_pair(i);
            break;
        }
    }
    
    itype v3_sin_pair= (v3_sin_pair_id!=-1) ? mesh.get_triangle(t1).TV(v3_sin_pair_id):-1;
    // if(v3_sin_pair==v2)
    //     {
    //        // cout<<"v3_sin:"<<v3_sin<<" v2:"<<v2<<endl;
    //       //  cout<<"v3 sin pair is v2"<<endl;  // why this case? 
    //   //      return false;
    //     }

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

    itype t3_des=-1;
    itype t3_sin=-1;
    bool v2_is_critical=true;
   
   // boost::dynamic_bitset<> edge_is_critical(vt2.size(),true);
    map<int,ivect> ets;
    for(int i=0; i<vt2.size(); i++){
        
        Triangle &t = mesh.get_triangle(vt2[i]);
        int v2_pos=t.vertex_index(v2);
        for(int j=0;j<t.vertices_num();j++){
            if(j!=v2_pos)
                {
                itype v_adj=t.TV(3-(v2_pos+j));
                ets[v_adj].push_back(vt2[i]);
                    //vv2.insert(t.TV(j));
                // t.TE(j,new_e);
                // v2_edges.insert(new_e);            
            //     if(!gradient.is_edge_critical(new_e,vt2[i],mesh)){
            //         edge_is_critical[v_adj] = 0;
            //   //  if(v2==217)
            //   //  cout<<"edge "<<new_e[0]<<", "<<new_e[1]<<" is not critical"<<endl;
            //     }
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

        if(gradient.is_triangle_critical(vt2[i])) {

        //     cout<<"vt2 is critical"<<endl;
            return false;}
        //Instead of searching for vtstar, we check all the triangles here
        if(!gradient.is_vertex_critical(v2,vt2[i],mesh)) 
        v2_is_critical=false;
    }
cout<<"breakpoint2"<<endl;
    if(v2_is_critical)
     {
       //  cout<<"v2 is critical"<<endl;
         return false;}

    for(auto it=ets.begin();it!=ets.end();it++){
        ivect e={it->first,v2};
        itype et1=it->second[0];
        itype et2=it->second[1];
    
        cout<<it->first<<": "<<et1<<", "<<et2<<endl;
        if(gradient.is_edge_critical(e,et1,mesh)&&gradient.is_edge_critical(e,et2,mesh))
        {
        //    cout<<"vv(v2) has critical edge"<<endl;
            return false;
        }
    }
cout<<"breakpoint1"<<endl;
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
            //cout<<"vt1 is critical"<<endl;
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
      //   cout<<"edge is critical"<<endl;
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

 //   cout<<"v1: "<<v1<<" v1_pair:"<<v1_pair<<endl;
  //  cout<<"v2: "<<v2<<" v2_pair:"<<v2_pair<<endl;
    if(v1_pair!=v2&&v2_pair!=v1){
  //      cout<<"edge is not paired with v1 or v2"<<endl;

        return false;
    }

    if(v3_sin_pair!=-1&& v3_sin_pair==v2){
        // cout<<"v3_sin is paired with v2"<<endl;
        // cout<<"t3_adj_sin:"<<t3_adj_sin<<endl;
     gradient.update_VE_adj_T(t3_adj_sin,v3_sin,v1,mesh,gradient);   
    }
    else if (v3_sin_pair!=-1&& v3_sin_pair==v1)
    {
        //    cout<<"v3_sin is paired with v1"<<endl;
        //    cout<<"v3_sin:"<<v3_sin<<" v3_sin_pair:"<<v3_sin_pair<<endl;
        //  cout<<"t3_sin:"<<t3_sin<<endl;
     gradient.update_VE_adj_T(t3_sin,v3_sin,v2,mesh,gradient);
    }
    else if (v1_pair!=-1&&v1_pair==v3_sin)
    {
        // cout<<"v1 is paired with v3_sin"<<endl;
        // cout<<"v1:"<<v1<<" v1_pair:"<<v1_pair<<endl;
     gradient.update_VE_adj_T(t3_sin,v2,v3_sin,mesh,gradient);
    }
    else if(v1_pair==v2&&v2_pair==v3_sin)
    {
        // cout<<"v2 is paired with v3_sin"<<endl;
        // cout<<"v2:"<<v2<<" v2_pair:"<<v2_pair<<endl;
        gradient.update_VE_adj_T(t3_adj_sin,v1,v3_sin,mesh,gradient);
    }

    if(v3_des_pair!=-1&& v3_des_pair==v2){

        // cout<<"v3_des is paired with v2"<<endl;
        //  cout<<"t3_adj_des:"<<t3_adj_des<<endl;
     gradient.update_VE_adj_T(t3_adj_des,v3_des,v1,mesh,gradient);   
    }
    else if (v3_des_pair!=-1&& v3_des_pair==v1)
    {
        // cout<<"v3_des is paired with v1"<<endl;
        //  cout<<"t3_des:"<<t3_des<<endl;
      //  gradient.update_VE_adj_T(t3,v3_des,v2,mesh);
     gradient.update_VE_adj_T(t3_des,v3_des,v2,mesh,gradient);
    }
    else if (v1_pair!=-1&&v1_pair==v3_des)
    {
        //  cout<<"v1 is paired with v3_des"<<endl;

        // cout<<"v1:"<<v1<<" v1_pair:"<<v1_pair<<endl;
     gradient.update_VE_adj_T(t3_des,v2,v3_des,mesh,gradient);
    }
    else if(v1_pair==v2&&v2_pair==v3_des)
    {
        // cout<<"v2 is paired with v3_des"<<endl;
        // cout<<"v2:"<<v2<<" v2_pair:"<<v2_pair<<endl;
        gradient.update_VE_adj_T(t3_adj_des,v1,v3_des,mesh,gradient);
    }
    return true;
    //Triangle 


}





void Gradient_Aware_Simplifier::get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1,bool& v1_is_border, bool& v2_is_border, Node_V *&outer_v_block,
                                                Node_V &n, Mesh &mesh, leaf_VT &vts,boost::dynamic_bitset<>is_border_edge, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params,PRT_Tree &tree){

    //cout<<"[NOTICE]get edge relation"<<endl;
    outer_v_block = NULL;
    /// inverted order as I only need the block indexing v1
    if(e[1]>e[0]){
    vt1 = Contraction_Simplifier::get_VT(e[1],n,mesh,vts,cache,tree,outer_v_block,params);
    vt0 = Contraction_Simplifier::get_VT(e[0],n,mesh,vts,cache,tree,outer_v_block,params);
   
    v2_is_border=is_border_edge[e[1]-n.get_v_start()];
    if(n.indexes_vertex(e[0])){
        v1_is_border=is_border_edge[e[0]-n.get_v_start()];
    }
    else{
        for(auto it=vt0->begin();it!=vt0->end();it++){
            Triangle& t=mesh.get_triangle(*it);
            int v_pos=t.vertex_index(e[0]);
            for(int v1=1;v1<t.vertices_num();v1++){
                if(t.is_border_edge((v1+v_pos)%t.vertices_num()))
                {
                    v1_is_border=true;
                    break;
                }
            }

        }
    }
   
   
    }
    else{

    vt0 = Contraction_Simplifier::get_VT(e[0],n,mesh,vts,cache,tree,outer_v_block,params);
    vt1 = Contraction_Simplifier::get_VT(e[1],n,mesh,vts,cache,tree,outer_v_block,params);

    v1_is_border=is_border_edge[e[0]-n.get_v_start()];
    if(n.indexes_vertex(e[1])){
        v2_is_border=is_border_edge[e[1]-n.get_v_start()];
    }
    else{
        for(auto it=vt1->begin();it!=vt1->end();it++){
            Triangle& t=mesh.get_triangle(*it);
            int v_pos=t.vertex_index(e[1]);
            for(int v1=1;v1<t.vertices_num();v1++){
                if(t.is_border_edge((v1+v_pos)%t.vertices_num()))
                {
                    v2_is_border=true;
                    break;
                }
            }

        }
    }


    }

    Contraction_Simplifier::get_ET(e,et,n,mesh,vts);

}