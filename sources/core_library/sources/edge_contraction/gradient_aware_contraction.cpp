#include "gradient_aware_contraction.h"

void Gradient_Aware_Simplifier::gradient_aware_simplify(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli,Forman_Gradient &gradient)   
{

    cerr<<"[NOTICED] Cache size: "<<cli.cache_size<<endl;
    LRU_Cache<int,leaf_VT> cache(cli.cache_size); // the key is v_start while the value are the VT relations
    contraction_parameters params;
    if(cli.contract_all_edges==true)
       params.contract_all_possible_edges();
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
time.start();
        trianglePlane =vector<dvect>(mesh.get_triangles_num(),dvect(4,0));
        initialQuadric = vector<Matrix>(mesh.get_vertices_num()+1,Matrix(0.0));
        cout<<"=========Calculate triangle plane========"<<endl;
        compute_triangle_plane(mesh,trianglePlane);
        cout<<"=========Calculate initial QEM========"<<endl;
        compute_initial_QEM(mesh,trianglePlane);
        vector<dvect>().swap(trianglePlane);
    time.stop();
    time.print_elapsed_time("[TIME] Calculating initial QEM: ");
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
    Gradient_Aware_Simplifier::update_mesh_and_tree(tree,mesh,params,gradient);
}

void Gradient_Aware_Simplifier::gradient_aware_simplify_parallel(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli,Forman_Gradient &gradient)   
{
    contraction_parameters params;
    if(cli.contract_all_edges==true)
       params.contract_all_possible_edges();
    params.set_maximum_limit(cli.maximum_limit);

    
    if(cli.QEM_based)
        params.queue_criterion_QEM();
    else
    {
        params.queue_criterion_length();
    }
    // Set to be parallel mode
    params.parallel_compute();
    Timer time;
    int simplification_round;
    int round = 1;

 time.start();
    cout<<"Number of threads used in the simplification:"<<omp_get_max_threads()<<endl;
  //  const int t_num = mesh.get_triangles_num();
    const int v_num = mesh.get_vertices_num();
    const int l_num = tree.get_leaves_number();
    // omp_lock_t lock[t_num];
  //  t_locks.resize(t_num);
    v_locks.resize(v_num);
    l_locks.resize(l_num);

// #pragma omp parallel for
//     for (int i = 0; i < v_num; i++)
//         omp_init_lock(&(v_locks[i]));
//     cout << "Initialize v_locks" << endl;

#pragma omp parallel for
    for (int i = 0; i < l_num; i++)
        omp_init_lock(&(l_locks[i]));
    cout << "Initialize l_locks" << endl;
    time.stop();
      time.print_elapsed_time("[TIME] Initialization of locks:  ");



    if(params.is_QEM()){

   //     trianglePlane =vector<dvect>(mesh.get_triangles_num(),dvect(4,0));
        initialQuadric = vector<Matrix>(mesh.get_vertices_num()+1,Matrix(0.0));
        time.start();
        cout<<"=========Calculate triangle plane========"<<endl;
    //    compute_triangle_plane(mesh,trianglePlane);
      //  time.stop();
        // time.print_elapsed_time("[TIME] Calculating Planes: ");
        // time.start();
        cout<<"=========Calculate initial QEM========"<<endl;
     //  compute_initial_QEM_parallel(tree,mesh,trianglePlane);
       compute_initial_plane_and_QEM_parallel(tree,mesh);
        time.stop();
        time.print_elapsed_time("[TIME] Calculating initial QEM: ");
      //  vector<dvect>().swap(trianglePlane);
        cerr << "[MEMORY] peak for computing QEM: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

        }
   time.start();
    while(1)
    {
        simplification_round = params.get_contracted_edges_num();  //checked edges

        /// HERE YOU NEED TO DEFINE A PROCEDURE FOR SIMPLIFY THE TIN BY USING THE SPATIAL INDEX
        this->simplify_compute_parallel(mesh, tree.get_subdivision(), params, tree,gradient);

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
   
        ///UPDATE: update data structure and conflict nodes lists after each round
    vector<int>().swap(v_in_leaf);
    lists_leafs().swap(conflict_leafs);

    
    
   // cout<<"number of remaining triangles: "<<tree.get_mesh().get_triangles_num()<<endl;

        // for (unsigned i = 0; i < tree.get_leaves_number(); i++)
        // {
        //     Node_V * leaf = tree.get_leaf(i);
        //                  ///TEMP
        //         for(auto it=leaf->get_t_array_begin();it!=leaf->get_t_array_end();it++)
        //         cout<<*it<<", ";
        //         cout<<endl;
        // }


        if(simplification_round == params.get_contracted_edges_num())
            break;
    Contraction_Simplifier::preprocess(tree,mesh,cli);
    }

// #pragma omp parallel for
//     for (int i = 0; i < v_num; i++)
//         omp_destroy_lock(&(v_locks[i]));

#pragma omp parallel for
    for (int i = 0; i < l_num; i++)
        omp_destroy_lock(&(l_locks[i]));

    time.stop();

    ///// Clear all the auxiliary data structures.
    //v_locks.clear();
    vector<omp_lock_t>().swap(v_locks);
    vector<omp_lock_t>().swap(l_locks);
    // vector<int>().swap(v_in_leaf);
    // lists_leafs().swap(conflict_leafs);
   // l_locks.clear();
    

    if(!cli.debug_mode)
        time.print_elapsed_time("[TIME] Edge contraction simplification: ");
    cerr << "[MEMORY] peak for Simplification: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

  //  cerr << "[RAM peak] for contracting a simplicial complex: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

    /// finally we have to update/compress the mesh and the tree
    // Gradient_Aware_Simplifier::update_mesh_and_tree(tree,mesh,params,gradient);
    // cerr << "[MEMORY] peak for mesh and tree updating: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
    time.start();
    Gradient_Aware_Simplifier::update_mesh_and_tree(tree,mesh,params,gradient);
    time.stop();
    time.print_elapsed_time("[TIME] Mesh and tree updating: ");

    cerr << "[MEMORY] peak for mesh and tree updating: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
 
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
    map<vector<int>, double> updated_edges;
    int edge_num=edges.size();
    int edges_contracted_leaf=0;
    //cout<<"Edge number:"<<edges.size()<<endl;
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
        //DISABLED GRADIENT CHECK FOR NOW TO CHECK THE CORRECTNESS OF PARALLEL COMPUTATION
        if(link_condition(e[0],e[1],*vt0,*vt1,et,mesh)/*&&valid_gradient_configuration(e[0],e[1],*vt0,*vt1,et,v1_is_border,v2_is_border,gradient,mesh)*/){
        contract_edge(e,et,*vt0,*vt1,*outer_v_block,edges,n,mesh,params,gradient,updated_edges);
        edges_contracted_leaf++;
    // break;
        }
    // cout<<"Number of edges remaining:"<<edges.size()<<endl;
    delete vt0,vt1,outer_v_block;
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
map<vector<int>, double> updated_edges;
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
        if(link_condition(e[0],e[1],*vt0,*vt1,et,mesh)/*&&valid_gradient_configuration(e[0],e[1],*vt0,*vt1,et,v1_is_border,v2_is_border,gradient,mesh)*/){
        contract_edge(e,et,*vt0,*vt1,*outer_v_block,edges,n,mesh,params,gradient,updated_edges);
        edges_contracted_leaf++;
    // break;
        }

    }

}

void Gradient_Aware_Simplifier::contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1,  Node_V &outer_v_block, edge_queue &edges,
                                           Node_V &n, Mesh &mesh, contraction_parameters &params,Forman_Gradient &gradient, map<vector<int>, double>& updated_edges)
{
 //   cout<<"[EDGE CONTRACTION] v1 and v2:"<<e[0]-1<<", "<<e[1]-1<<endl;
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
    //Contraction_Simplifier::update(e,vt0,vt1,n,outer_v_block,edges,mesh,params);

    if (!params.is_parallel())
        Contraction_Simplifier::update(e, vt0, vt1, n, outer_v_block, edges, mesh, params,updated_edges);
    else
    {

        Contraction_Simplifier::update_parallel(e, vt0, vt1, n, outer_v_block, edges, mesh, params,updated_edges);
    }

    // we remove v2 and the triangles in et
    Contraction_Simplifier::remove_from_mesh(e[1],et,mesh,params); 
    // finally we clear the VT(v2)
    vt1.clear();
    //et.clear();
}


bool Gradient_Aware_Simplifier::valid_gradient_configuration(int v1,int v2, VT &vt1, VT &vt2,ET& et ,bool v1_is_border, bool v2_is_border, Forman_Gradient &gradient, Mesh &mesh){

    bool debug=false;

    // debug=true;
    if(v1==23423&&v2==23422)
        debug =true;
    if(v1==23415&&v2==23517)
        debug =true;

  if(debug)  
    cout<<"[debug]checking edge "<<v1<<", "<<v2<<endl;
    if(v1_is_border||v2_is_border){
        if(debug)
      cout<<"border edge"<<endl;
    return false;}
    if(vt1.size()<4||vt2.size()<4){
       if(debug)
        cout<<"less than 4 triangles"<<endl;
        return false;
    }
    
    int t1=et.first;
    int t2=et.second;
    if(gradient.is_triangle_critical(t1)||gradient.is_triangle_critical(t2))
    {
        if(debug)
        cout<<"t1 or t2 is critical"<<endl;
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
            cout<<*it<<endl;
            cout<<t<<endl;

        }
        cout<<endl;

        cout<<"[DEBUG]vt2 v2:"<<v2<<endl;
        for(auto it=vt2.begin();it!=vt2.end();it++){
             Triangle t=mesh.get_triangle(*it);
             cout<<*it<<endl;
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
        if(debug)
            cout<<"vt2 is critical"<<endl;
            return false;}
        //Instead of searching for vtstar, we check all the triangles here
        if(!gradient.is_vertex_critical(v2,vt2[i],mesh)) 
        v2_is_critical=false;
    }

    if(v2_is_critical)
     {
         if(debug)
        cout<<"v2 is critical"<<endl;
         return false;}

    for(auto it=ets.begin();it!=ets.end();it++){
        ivect e={it->first,v2};
        if(debug)
        cout<<"ets size"<<it->second.size()<<endl;
        itype et1=it->second[0];
        itype et2=it->second[1];
    
      //  cout<<it->first<<": "<<et1<<", "<<et2<<endl;
        if(gradient.is_edge_critical(e,et1,mesh)&&gradient.is_edge_critical(e,et2,mesh))
        {
            if(debug){
            cout<<"ET: "<<et1<<", "<<et2<<endl;
            cout<<mesh.get_triangle(et1)<<endl;
            cout<<mesh.get_triangle(et2)<<endl;
            cout<<"vv(v2) has critical edge"<<endl;
            }
            return false;
        }
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
            if(debug)
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
         if(debug)
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

 //   cout<<"v1: "<<v1<<" v1_pair:"<<v1_pair<<endl;
  //  cout<<"v2: "<<v2<<" v2_pair:"<<v2_pair<<endl;
    if(v1_pair!=v2&&v2_pair!=v1){
        if(debug)
       cout<<"edge is not paired with v1 or v2"<<endl;

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

       if(debug)
    cout<<"valid gradient condition"<<endl;
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



void Gradient_Aware_Simplifier::simplify_compute_parallel(Mesh &mesh,  Spatial_Subdivision &division, contraction_parameters &params, PRT_Tree &tree, Forman_Gradient &gradient)
{

    //boost::dynamic_bitset<> conflict_nodes(tree.get_leaves_number());
    ivect nodes_status(tree.get_leaves_number(), 0);
    int processed_node = 0;
    bool processed = false;
    do
    {
       // processed = false;
        
#pragma omp parallel for reduction(+:processed_node)// schedule(dynamic,1)
        for (unsigned i = 0; i < tree.get_leaves_number(); i++)
        {
            Node_V *leaf = tree.get_leaf(i);
            if(nodes_status[i]==-1)
               continue;
            //cout<<*leaf<<endl;
        if (!leaf->indexes_vertices())
         {
           processed_node++;
           nodes_status[i] =-1;
            continue;
          }
            //check the array of conflict_nodes
            // if nodes_status[i]==1, then continue
          //  cout << "Current leaf node:" << i << " on thread " << omp_get_thread_num() << endl;
    omp_set_lock(&(l_locks[i]));
            if (nodes_status[i] != 0)
            {
                omp_unset_lock(&(l_locks[i]));
                //   cout<<"Node "<<i<<" is skipped."<<endl;
                continue;
            }
            // else we set the conflict nodes of leaf[i] to be 1 in nodes_status
            else
            {
                nodes_status[i]=2;
                omp_unset_lock(&(l_locks[i]));
                iset conflicts = conflict_leafs[i];
                // Check if the conflict nodes were set to 1 already
               // bool shared_conflicts = false;
               bool cannot_process = false;
                for (iset_iter it = conflicts.begin(); it != conflicts.end(); it++)
                {
                    //  cout<<"set leaf node:"<<*it<<" on thread "<<omp_get_thread_num()<<endl;
                    omp_set_lock(&(l_locks[*it])); // leafs should be locked when being checked and updated.
                    int status = 0;
                    //   #pragma omp atomic read
                    status = nodes_status[*it];
                    if (status == 1||status==2) //it is conflicted with another node being processed
                    {
                        // cout<<"conflict node id:"<<*it<<" with "<<i<<" on thread "<<omp_get_thread_num()<<endl;
                        omp_unset_lock(&(l_locks[*it]));
                        //it++; No need, the current one should not be change to 0
                        for (iset_iter it2 = conflicts.begin(); it2 != it; it2++)
                        {
                            //     cout<<"unset leaf node:"<<*it2<<" on thread "<<omp_get_thread_num()<<endl;
                            omp_set_lock(&(l_locks[*it2]));
                            if (nodes_status[*it2] == 1)
                                nodes_status[*it2] = 0;
                            omp_unset_lock(&(l_locks[*it2]));
                        }
                        //omp_unset_lock(&(l_locks[*it]));
                        // unset the locks that have been set
                        cannot_process = true;
                        //cout<<"neighbor conflict"<<endl;

                        break;
                    }
                    else if (status == 0)
                    {
                        nodes_status[*it] = 1;
                    }
                    omp_unset_lock(&(l_locks[*it]));
                }
                if (cannot_process == true)
                {
                    omp_set_lock(&(l_locks[i]));
                    nodes_status[i]=0;
                    omp_unset_lock(&(l_locks[i]));
                    continue;
                }
  
                
                //  cout<<"Node "<<i<<" will be processed."<<endl;
                // omp_set_lock(&(l_locks[i]));
                // // set nodes_status[i]=2 when node is being processed
                // nodes_status[i] = 2;
                // omp_unset_lock(&(l_locks[i]));
                processed_node=processed_node+1;
              //  processed = true;
              // cout << "Start simplification" << endl;
               if(params.is_QEM()==true)
                simplify_leaf_cross_QEM(*leaf, i, mesh, params, tree,gradient);
               else
                simplify_leaf_cross(*leaf, i, mesh, params, tree,gradient);

              //  cout << "Finish simplification" << endl;
                //set nodes_status[i]=-1 after processing

                omp_set_lock(&(l_locks[i]));
                nodes_status[i] = -1;
                omp_unset_lock(&(l_locks[i]));
                //cout<<"unset leaf node:"<<i<<" on thread "<<omp_get_thread_num()<<endl;

                for (iset_iter it = conflicts.begin(); it != conflicts.end(); it++)
                {
                    //  cout<<" set leaf lock "<<*it<<" on thread "<<omp_get_thread_num()<<endl;
                    omp_set_lock(&(l_locks[*it]));
                    int status = 0;
                    //   #pragma omp atomic read
                    status = nodes_status[*it];
                    if (status == 1)
                    {
                        nodes_status[*it] = 0;
                    }
                    //cout<<"unset leaf node:"<<*it<<" on thread "<<omp_get_thread_num()<<endl;
                    omp_unset_lock(&(l_locks[*it]));
                }
            }
        }
       // cout << "finished one for loop" << endl;
 //cerr << "[MEMORY] peak for a simplification round:" << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
//cerr<<"Number of processed nodes:"<<processed_node<<endl;
    //} while (processed == true);
    } while (processed_node!=tree.get_leaves_number());
}


void Gradient_Aware_Simplifier::simplify_leaf_cross(Node_V &n, int n_id, Mesh &mesh, contraction_parameters &params, PRT_Tree &tree,Forman_Gradient &gradient)
{

    if (!n.indexes_vertices())
        return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();
    itype v_range = v_end - v_start; 

    //cout<<"Simplification in leaf."<<endl;
    // leaf_VT local_vts(v_range,VT());
    // n.get_VT(local_vts,mesh);
    boost::dynamic_bitset<> is_v_border(v_end-v_start);
    leaf_VT local_vts(v_range, VT());
    n.get_VT_and_border(local_vts,is_v_border, mesh);

    // Create a priority queue of candidate edges
    edge_queue edges;
    // if (params.is_parallel())
    // {
    //     find_candidate_edges_parallel(n, mesh, local_vts, edges, params, true);
    // }
    // else
    
    find_candidate_edges(n, mesh, local_vts, edges, params);
map<vector<int>, double> updated_edges;
    int edge_num = edges.size();
    int edges_contracted_leaf = 0;
    
    // cout<<"Number of threads used: "<<omp_get_num_threads()<<endl;
    // cout<<"Current thread id: "<<omp_get_thread_num()<<endl;

    params.add_edge_queue_size(edges.size());
    map<int, leaf_VT> local_cache;
    while (!edges.empty())
    {
        Geom_Edge *current = edges.top();
        ivect e = current->edge;
            //cout<<"Start contraction."<<endl;
        //  cout<<"Edge Length:"<<current->val<<endl;

        edges.pop();

        if (mesh.is_vertex_removed(e[0]) || mesh.is_vertex_removed(e[1]))
        {

            delete current;
            continue;
        }

        ET et(-1, -1);
        VT *vt0 = NULL, *vt1 = NULL;
        Node_V *outer_v_block = NULL;
        bool v1_is_border=false, v2_is_border=false;
      
        Contraction_Simplifier::get_edge_relations(e, et, vt0, vt1,v1_is_border,v2_is_border, outer_v_block, n, mesh, local_vts, is_v_border, local_cache, params, tree);

        // if(params.is_parallel()){
        VV vv_locks;
        if (link_condition(e[0], e[1], *vt0, *vt1, et, n, *outer_v_block, vv_locks, mesh)&&not_fold_over(e[0], e[1], *vt0, *vt1, et, mesh)&&valid_gradient_configuration(e[0],e[1],*vt0,*vt1,et,v1_is_border,v2_is_border,gradient,mesh))
        {
            contract_edge(e, et, *vt0, *vt1, *outer_v_block, edges, n, mesh, params,gradient,updated_edges);
            edges_contracted_leaf++;
            // break;

            // A new step for cross edge case
            // Check possible new conflict nodes by checking the vv_locks
            // vv_locks stores all the vertices in the VV(v0) & VV(v1) that are not contained by n or outer_v_block
            update_conflict_nodes(vv_locks, n_id, tree);
        }
        // for (iset_iter it = vv_locks.begin(); it != vv_locks.end(); it++)
        // {
        //     omp_unset_lock(&(v_locks[*it - 1]));
        // }
        // }
        // else
        // {
        // if (link_condition(e[0], e[1], *vt0, *vt1, et, mesh))
        // {
        //     contract_edge(e, et, *vt0, *vt1, *outer_v_block, edges, n, mesh, cache, params);
        //     edges_contracted_leaf++;
        //     // break;
        // }
        // }
        delete current;
    }

    // leaf_VV vvs;
    // n.get_VV(vvs,mesh);

}

void Gradient_Aware_Simplifier::simplify_leaf_cross_QEM(Node_V &n, int n_id, Mesh &mesh, contraction_parameters &params, PRT_Tree &tree,Forman_Gradient &gradient)
{

    if (!n.indexes_vertices())
        return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();
    itype v_range = v_end - v_start;

    //cout<<"Simplification in leaf."<<endl;
    boost::dynamic_bitset<> is_v_border(v_end-v_start);
    leaf_VT local_vts(v_range, VT());
    n.get_VT_and_border(local_vts,is_v_border, mesh);



    // Create a priority queue of candidate edges
    edge_queue edges;

    find_candidate_edges_QEM(n, mesh, local_vts, edges, params);
    map<vector<int>, double> updated_edges;
    int edge_num = edges.size();
    int edges_contracted_leaf = 0;
   // cout << "Edge number:" << edges.size() << endl;
    // cout<<"Number of threads used: "<<omp_get_num_threads()<<endl;
    // cout<<"Current thread id: "<<omp_get_thread_num()<<endl;

    params.add_edge_queue_size(edges.size());
    map<int, leaf_VT> local_cache;
    while (!edges.empty())
    {
        Geom_Edge *current = edges.top();
        ivect e = current->edge;
        //    cout<<"Start contraction."<<endl;
     //   cout<<"Edge Length:"<<current->val<<endl;

        edges.pop();

        if (mesh.is_vertex_removed(e[0]) || mesh.is_vertex_removed(e[1]))
        {

            delete current;
            continue;
        }
        ivect sorted_e = e;
        sort(sorted_e.begin(), sorted_e.end());
        auto it = updated_edges.find(sorted_e);
        if (it != updated_edges.end())
        {
            //int tmp=-1;
            // double error = compute_error(e[0],e[1],mesh,tmp);
            if (fabs(it->second - current->val) > Zero)
            {
                delete current;
                continue;
            }
        }
        ET et(-1, -1);
        VT *vt0 = NULL, *vt1 = NULL;
        Node_V *outer_v_block = NULL;
        bool v1_is_border=false, v2_is_border=false;

        Contraction_Simplifier::get_edge_relations(e, et, vt0, vt1,v1_is_border,v2_is_border, outer_v_block, n, mesh, local_vts, is_v_border, local_cache, params, tree);
             if(e[0]==23423&&e[1]==23422)
        {cout<<"Error value of 23422, 23421:"<<current->val<<endl;
        printf("%.10lf ", current->val);
        }
        if(e[0]==23421&&e[1]==23422)
        {cout<<"Error value of 23420, 23421:"<<endl;
        printf("%.10lf ", current->val);
        }        // if(params.is_parallel()){
        VV vv_locks;
        if (link_condition(e[0], e[1], *vt0, *vt1, et, n, *outer_v_block, vv_locks, mesh)&&not_fold_over(e[0], e[1], *vt0, *vt1, et, mesh)&&valid_gradient_configuration(e[0],e[1],*vt0,*vt1,et,v1_is_border,v2_is_border,gradient,mesh))
        {
            contract_edge(e, et, *vt0, *vt1, *outer_v_block, edges, n, mesh, params,gradient,updated_edges);
            edges_contracted_leaf++;
            // break;

             n_id=v_in_leaf[e[0]];
            // A new step for cross edge case
            // Check possible new conflict nodes by checking the vv_locks
            // vv_locks stores all the vertices in the VV(v0) & VV(v1) that are not contained by n or outer_v_block
            update_conflict_nodes(vv_locks, n_id, tree);
        }
        // for (iset_iter it = vv_locks.begin(); it != vv_locks.end(); it++)
        // {
        //     omp_unset_lock(&(v_locks[*it - 1]));
        // }
    delete current;
   // delete vt0,vt1,outer_v_block;
    }

    // leaf_VV vvs;
    // n.get_VV(vvs,mesh);
}


void Gradient_Aware_Simplifier::update_mesh_and_tree(PRT_Tree &tree, Mesh &mesh, contraction_parameters &params, Forman_Gradient &gradient)
{
    // Timer time;

    ///  UPDATE OF MESH AND TREE
    ivect new_v_positions;
    ivect new_t_positions;
    ivect surviving_vertices;

   // time.start();
    //    cerr<<"[TREE] compact vertices lists"<<endl;
    tree.compact_vertices_lists(tree.get_root(), mesh, surviving_vertices);
    // time.stop();
    // time.print_elapsed_time("[TIME] Compact tree vertices lists: ");
    // cerr << "[MEMORY] peak for compacting tree vertices lists: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    //    print_container_content("surviving vertices: ",surviving_vertices);
    //    mesh.print_mesh(cout);
    //    int a; cin>>a;

    // time.start();
    //    cerr<<"[MESH] compact"<<endl;
    Mesh_Updater mu;
    cout << "number of surviving vertices:" << surviving_vertices.size() << endl;
    mu.clean_vertices_array(mesh, new_v_positions, surviving_vertices);
    //cout << "number of deleted triangles:" << params.get_counter() << endl;

    gradient.reorder_forman_gradient(mesh);
    // if(params.is_QEM()){
    //     update_QEM(surviving_vertices,mesh);
    // }

    /// NEW: the update_and_compact procedure check internally if we have removed all the top d-simplices
    bool all_deleted = mu.update_and_clean_triangles_arrays(mesh, new_v_positions, new_t_positions, params.get_counter());
    // time.stop();
    // time.print_elapsed_time("[TIME] Compact and update mesh: ");
    // cerr << "[MEMORY] peak for compacting and updating the mesh: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;



    // cerr<<"[STAT] mesh "<<endl;
    // cerr<<"------------"<<endl;
    // mesh.print_mesh_stats(cerr);
    // cerr<<"------------"<<endl;

    // time.start();
    //    cerr<<"[TREE] update indices in the tree"<<endl;
    ///TODO: Check triangle intersection before updating the tree.

    tree.update_tree(tree.get_root(), new_v_positions, new_t_positions, all_deleted,1);
    // time.stop();
    // time.print_elapsed_time("[TIME] Update tree (top-simplices): ");
    // cerr << "[MEMORY] peak for updating the tree (top-simplices): " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    //    Reindexer r;
    //    r.reorganize_index_and_mesh(tree,mesh,false);
// for(int i = 0; i<tree.get_leaves_number();i++)
// {
//     Node_V * l = tree.get_leaf(i);
//     cout<<"Leaf "<<i<<endl;
//     cout<<*l<<endl;
//     for(auto it = l->v_array_begin_iterator();it!=l->v_array_end_iterator();it++)
//     cout<<*it<<", ";
//     cout<<endl;

// }



    //cerr << "[RAM peak] for updating the mesh and the tree: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
}