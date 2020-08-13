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

#include "forman_gradient_simplifier.h"

void Forman_Gradient_Simplifier::exec_global_topological_simplification(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,
                                                                        OpType operation, int cache_size, coord_type persistence)
{
    /// reset the statistical counters
    reset_extraction_critica_counters();

    mig_cache cache = mig_cache(cache_size);
    ig_paths paths;

    Timer time;
    time.start();

    forman_ig.init(); /// init again the MIG structures
    this->incidence_graph_extraction(n,mesh,gradient,division,n,OUTPUT,cache,paths); /// we need to encode explicitly the MIG
    if(!paths.visited_all())
    {
        cout<<"[ERROR] the labels are not correctly assigned."<<endl;
    }
    cache.reset();
    print_arc_nums();
    this->global_topological_simplification(n,mesh,gradient,division,n,operation,cache,persistence);
    time.stop();
    time.print_elapsed_time("TOT MIG + global simplification time: ");

    /// after the global simplification we check if we have some MIG arc that have a persistence below the target
    /// and could be simplified (i.e. saddle arcs == 2)
    /// for debug
    int not_simpl = 0;
    for(int i=0; i<2; i++)
    {
        for(set<Arc*>::iterator it = forman_ig.getLevelArcs(i).begin(); it != forman_ig.getLevelArcs(i).end(); ++it)
        {
            if((*it)->getLabel() == 1)
            {
                Vertex &v1 = mesh.get_vertex((*it)->getNode_i()->get_critical_index());

                if(i==1)
                {
                    // ivect t;
                    // mesh.get_triangle((*it)->getNode_j()->get_critical_index()).convert_to_vec(t);
                    
                    Vertex &v2 = mesh.get_vertex(get_max_elevation_vertex(mesh.get_triangle((*it)->getNode_j()->get_critical_index())));
                    if(abs(v1.get_z() - v2.get_z()) <= persistence)
                    {
                        iNode *saddle = (iNode*)(*it)->getNode_i();
                        if(saddle->getArcs(false).size() == 2)
                        {
                            //cout<<i<<") pers["<<fabs(v1.get_z() - v2.get_z())<<"] "<<**it<<endl;
                            not_simpl++;
                        }
                    }
                }
                else{
                    Vertex &v2 = mesh.get_vertex((*it)->getNode_j()->get_critical_index());
                    if(abs(v1.get_z() - v2.get_z()) <= persistence)
                    {
                        iNode *saddle = (iNode*)(*it)->getNode_j();
                        if(saddle->getArcs(true).size() == 2)
                        {
                           // cout<<i<<") pers["<<fabs(v1.get_z() - v2.get_z())<<"] "<<**it<<endl;
                            not_simpl++;
                        }
                    }
                }

            }
        }
    }
    if(not_simpl > 0)
        cerr<<"GLOBAL SIMPLIFICATION WARNING: "<<not_simpl<<" arcs can be simplified, but are not.."<<endl;
}

void Forman_Gradient_Simplifier::global_topological_simplification(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,
                                                                   Node_V &root, OpType operation, mig_cache &cache, coord_type persistence)
{
    if (n.is_leaf())
    {
        /// if there are no vertices in the leaf we have nothing to do..
        if(!n.indexes_vertices())
            return;

        Timer time;
        local_VTstar_ET local_rels;

        if(operation == TIME_VERBOSE)
            time.start();
        Forman_Gradient_Topological_Relations::get_VTstar_ET(local_rels,n,mesh,gradient);
        if(operation == TIME_VERBOSE)
        {
            time.stop();
            rels_gathering += time.get_elapsed_time();
            time.start();
        }

//        map<pair<itype,itype>, iNode*> selle = forman_ig.getSaddles();
//        for(auto s : selle)
//            if(n.indexes_vertex(s.second->get_critical_index()))
//                cout<<"S("<<s.first.first<<" "<<s.first.second<<") -> "<<*(s.second)<<endl;

        this->global_topological_simplification_leaf(n,mesh,gradient,local_rels,cache,root,division,persistence);
        if(operation == TIME_VERBOSE)
        {
            time.stop();
            u_paths += time.get_elapsed_time();
        }

        if(operation == TIME_VERBOSE)
            time.start();
        //aggiungo alla cache la struttura della foglia corrente
        if(mesh.get_vertices_num() > n.get_v_end())
        {
            cache.addItem(n.get_v_end()+n.get_v_start(),local_rels);
        }
        if(operation == TIME_VERBOSE)
        {
            time.stop();
            cache_t += time.get_elapsed_time();
        }
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->global_topological_simplification(*n.get_son(i),mesh,gradient,division,root,operation,cache,persistence);
            }
        }
    }
}

void Forman_Gradient_Simplifier::global_topological_simplification_leaf(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, local_VTstar_ET &local_rels, mig_cache &cache,
                                                                  Node_V& root, Spatial_Subdivision &division, coord_type persistence)
{
//    cout<<n<<endl;
    priority_arcs_queue queue;
    this->build_persistence_queue(queue,local_rels.get_ETs(),mesh,gradient,persistence); /// add a subset of the arcs of the local MIG
    cout<<"======FINISHED BUILD QUEUE====="<<endl;
    if(max_priority_queue_size < (int)queue.size())
        max_priority_queue_size = queue.size();

    simplify(queue,n,mesh,gradient,local_rels, cache,root,division,persistence);
}

void Forman_Gradient_Simplifier::simplify(priority_arcs_queue &queue, Node_V &n, Mesh &mesh, Forman_Gradient &gradient, local_VTstar_ET &local_rels,
                                    mig_cache &cache, Node_V &root, Spatial_Subdivision &division, coord_type persistence)
{
    Topo_Sempl sempl;
    iNode* saddle=NULL;
    nNode* extrema=NULL;
    int lvl0=0,lvl1=0,lvl2=0;
    while(!queue.empty())
    {
        sempl = queue.top();
        queue.pop();

        if(sempl.arc->getLabel() != 1 )
        {
                        lvl2++;
            if(sempl.arc->getLabel() == -1)
                delete sempl.arc;
            continue;
        }

        if(sempl.lvl == 0)
        {
                        lvl0++;
            saddle  = (iNode*)sempl.arc->getNode_j();
            extrema = (nNode*)sempl.arc->getNode_i();
          //  cout<<"[Contraction] persistence value:"<<sempl.val<<" edge filtration: "<<sempl.filt_s0<<", "<<sempl.filt_s1<<"; extreme filtration: "<<sempl.filt_ex<<endl;
            if(saddle->getArcs(true).size() != 2) continue;
            // cout<<sempl.arc->getSimplexi()"Persistence value:"<<sempl.val<<endl;
            contraction(extrema, saddle, queue, mesh, gradient, local_rels, cache, n, root, division, persistence);

            refined_topo++;
        }
        else
        {
            lvl1++;
            saddle  = (iNode*)sempl.arc->getNode_i();
            extrema = (nNode*)sempl.arc->getNode_j();

            if(saddle->getArcs(false).size() != 2) continue;
           // cout<<"[Removal] persistence value:"<<sempl.val<<" edge filtration: "<<sempl.filt_s0<<", "<<sempl.filt_s1<<endl;
           // cout<<"extreme filtration: "<<sempl.filt_ex[0]<<"; "<<sempl.filt_ex[1]<<"; "<<sempl.filt_ex[2]<<endl;

            removal(extrema,saddle, queue, mesh, gradient, local_rels, cache, n, root, division, persistence);

            refined_topo++;
        }

        if(max_priority_queue_size < (int)queue.size())
            max_priority_queue_size = queue.size();
    }
   // cout<<"Level 0:"<<lvl0<<", Level 1:"<<lvl1<<", Level 2:"<<lvl2<<endl;
}

void Forman_Gradient_Simplifier::contraction(nNode *extrema, iNode *saddle, priority_arcs_queue &q, Mesh &mesh, Forman_Gradient &gradient, local_VTstar_ET &local_rels,
                                       mig_cache &cache, Node_V &n, Node_V &root, Spatial_Subdivision &division, coord_type persistence)
{
    vector<Arc*> arcs = saddle->get_vector_Arcs(true);

    pair<itype,itype> critical_edge_tetra=saddle->get_edge_id();
    ivect critical_edge;

    Triangle &first = mesh.get_triangle(critical_edge_tetra.first);

    if(critical_edge_tetra.second<0)
    {
        first.TE(-critical_edge_tetra.second-1,critical_edge);
    }
    else
    {
        Triangle &second = mesh.get_triangle(critical_edge_tetra.second);
        for(int i=0; i<3; i++)
        {
            if(!(second.has_vertex(first.TV(i))))
            {
                first.TE(i,critical_edge);
                break;
            }
        }
    }

    nNode* other_extrema;
    itype vertex,next_vertex;
    itype ending_path_simplex;
    //SETUp
    if(arcs[0]->getNode_i() == extrema){
        other_extrema = (nNode*)arcs[1]->getNode_i();
        next_vertex = first.TV(arcs[0]->getSimplexj());
        vertex = first.TV(arcs[1]->getSimplexj());
        ending_path_simplex = arcs[1]->getSimplexi();
    }
    else{
        other_extrema = (nNode*)arcs[0]->getNode_i();
        next_vertex = first.TV(arcs[1]->getSimplexj());
        vertex = first.TV(arcs[0]->getSimplexj());
        ending_path_simplex = arcs[0]->getSimplexi();
    }

    /// ----------------- GRADIENT UPDATES ----------------- ///
    contraction_update_gradient(vertex, extrema->get_critical_index(), next_vertex, critical_edge, mesh, gradient, local_rels, cache, n, root, division);

    /// ----------------- MIG UPDATES ----------------- ///
    arcs[0]->setLabel(-1);
    arcs[1]->setLabel(-1);
    /// remove maxima-saddle arcs
    remove_saddle_arcs(saddle,false,forman_ig);
    /// remove extrema arcs
    remove_extreme_arcs(extrema, saddle, other_extrema, ending_path_simplex, true, q, n, mesh, persistence);
    /// remove minima-saddle arcs
    remove_saddle_arcs(saddle,true,forman_ig);
    /// remove the saddle and the minimum
    // if((filtration[saddle->get_critical_index()-1]==260261)||(filtration[saddle->get_critical_index()-1]==260290)){
    // cout<<"CONTRACTION SADDLE: "<<*saddle<<endl;
    // cout<<"First triangle: ["<<mesh.get_triangle(saddle->get_edge_id().first).TV(0)<<", "<<mesh.get_triangle(saddle->get_edge_id().first).TV(1)<<", "<<mesh.get_triangle(saddle->get_edge_id().first).TV(2)<<endl;
    // cout<<"Second triangle: ["<<mesh.get_triangle(saddle->get_edge_id().second).TV(0)<<", "<<mesh.get_triangle(saddle->get_edge_id().second).TV(1)<<", "<<mesh.get_triangle(saddle->get_edge_id().second).TV(2)<<endl;
    // cout<<"DEBUG:"<<filtration[676950]<<"; "<<filtration[676951]<<"; "<<filtration[678751]<<";"<<filtration[677851]<<endl;

    // }
    cout<<"[CONTRACTION]Filtration value: "<<filtration[saddle->get_critical_index()-1]<<endl;
   // cout<<"[CONTRACTION]Minimum Filtration value: "<<filtration[extrema->get_critical_index()-1]<<endl;
    forman_ig.remove_saddle(saddle->get_edge_id(),saddle);
    forman_ig.remove_minimum(extrema->get_critical_index(),extrema);

//    ivect edge = {11376, 11377};
//    ET et = {22741, 22728};
//    cout<<"IS THE TARGET EDGE CRITICAL? "<<gradient.is_edge_critical(edge,et,mesh)<<endl;
}

void Forman_Gradient_Simplifier::removal(nNode *extrema, iNode *saddle, priority_arcs_queue &q, Mesh &mesh, Forman_Gradient &gradient, local_VTstar_ET &local_rels,
                                   mig_cache &cache, Node_V &n, Node_V &root, Spatial_Subdivision &division, coord_type persistence)
{
    vector<Arc*> arcs = saddle->get_vector_Arcs(false);
    nNode* other_extrema;
    itype triangle;
    itype ending_path_simplex;

    if(arcs[0]->getNode_j() == extrema){
        other_extrema = (nNode*)arcs[1]->getNode_j();
        triangle = arcs[0]->getSimplexi();
        ending_path_simplex = arcs[1]->getSimplexj();
    }
    else{
        other_extrema = (nNode*)arcs[0]->getNode_j();
        triangle = arcs[1]->getSimplexi();
        ending_path_simplex = arcs[0]->getSimplexj();
    }

    /// ----------------- GRADIENTE UPDATES ----------------- ///
    removal_update_gradient(triangle, saddle, mesh, gradient, local_rels.get_ETs(), cache, n, root, division);

    /// ----------------- MIG UPDATES ----------------- ///
    arcs[0]->setLabel(-1);
    arcs[1]->setLabel(-1);

    /// remove minima-saddle arcs
    remove_saddle_arcs(saddle,true,forman_ig);
    /// remove extrema arcs
    remove_extreme_arcs(extrema, saddle, other_extrema, ending_path_simplex, false, q, n, mesh, persistence);
    /// remove maxima-saddle arcs
    remove_saddle_arcs(saddle,false,forman_ig);
    /// remove the saddle and the maximum
   //cout<<"REMOVAL SADDLE: "<<*saddle<<endl;
    cout<<"[REMOVAL]Filtration value: "<<filtration[saddle->get_critical_index()-1]<<endl;
    //cout<<"[REMOVAL]Maximum Filtration value: "<<filtration[get_max_elevation_vertex(mesh.get_triangle(extrema->get_critical_index()))-1]<<endl;
    forman_ig.remove_saddle(saddle->get_edge_id(),saddle);
    forman_ig.remove_maximum(extrema->get_critical_index(),extrema);

//    ivect edge = {11376, 11377};
//    ET et = {22741, 22728};
//    cout<<"IS THE TARGET EDGE CRITICAL? "<<gradient.is_edge_critical(edge,et,mesh)<<endl;
}

void Forman_Gradient_Simplifier::remove_extreme_arcs(nNode* extrema, iNode* saddle, nNode *other_extrema, itype ending_path_simplex,
                                               bool is_minimum, priority_arcs_queue &q, Node_V &n, Mesh &mesh, coord_type persistence)
{
    iNode* node_saddle1=NULL;
    itype starting_path_simplex;
    coord_type val;

    set<Arc*>::iterator it=extrema->begin();
    while(it!=extrema->end())
    {
        if(is_minimum)
        {
            node_saddle1 = ((iNode*)(*it)->getNode_j());
            starting_path_simplex = (*it)->getSimplexj();
        }
        else
        {
            node_saddle1 = ((iNode*)(*it)->getNode_i());
            starting_path_simplex = (*it)->getSimplexi();
        }

        (*it)->setLabel(-1);

        forman_ig.removeArc(!is_minimum,*it);
        node_saddle1->removeArc(is_minimum, *it);
        extrema->removeArc(it);

        if(node_saddle1 != saddle)
        {
            Arc* existing_arc = forman_ig.already_connected(other_extrema,node_saddle1);
            if(existing_arc==NULL)
            {
                Arc* arco = NULL;

                if(is_minimum)
                    arco = forman_ig.addArc(other_extrema, ending_path_simplex, node_saddle1, starting_path_simplex, 0);
                else
                    arco = forman_ig.addArc(node_saddle1, starting_path_simplex, other_extrema, ending_path_simplex, 1);

//                cout<<"addArc: "<<*arco<<endl;

                /// the arc can be simplified AND the saddle is indexed into a leaf node that has already been visited
                /// if the saddle will be processed in a successive leaf node
                if(arco->getLabel() == 1 && n.visited_vertex(node_saddle1->get_critical_index()))
                {
                    itype index_i,index_j;
                    itype filt_i,filt_j;
                    ivect filt_ex;
                    /// check how to compute val..
                    if(is_minimum)
                    {
                        index_i=arco->getNode_i()->get_critical_index();
                        index_j=arco->getNode_j()->get_critical_index();
                        filt_ex.push_back(filtration[index_i-1]);
                        pair<itype,itype> critical_edge_tetra=((iNode*)arco->getNode_j())->get_edge_id();
                        ivect critical_edge;

                        Triangle &first = mesh.get_triangle(critical_edge_tetra.first);

                        if(critical_edge_tetra.second<0)
                        {
                            first.TE(-critical_edge_tetra.second-1,critical_edge);
                        }
                        else
                        {
                            Triangle &second = mesh.get_triangle(critical_edge_tetra.second);
                            for(int i=0; i<3; i++)
                            {
                                if(!(second.has_vertex(first.TV(i))))
                                {
                                    first.TE(i,critical_edge);
                                    break;
                                }
                            }
                        }
                         filt_i=(filtration[critical_edge[0]-1]>filtration[critical_edge[1]-1])?filtration[critical_edge[0]-1]:filtration[critical_edge[1]-1];
                        filt_j=(filtration[critical_edge[0]-1]>filtration[critical_edge[1]-1])?filtration[critical_edge[1]-1]:filtration[critical_edge[0]-1];

                    }

                    else
                        {index_i=arco->getNode_i()->get_critical_index();

                        pair<itype,itype> critical_edge_tetra=((iNode*)arco->getNode_i())->get_edge_id();
                        ivect critical_edge;

                        Triangle &first = mesh.get_triangle(critical_edge_tetra.first);

                        if(critical_edge_tetra.second<0)
                        {
                            first.TE(-critical_edge_tetra.second-1,critical_edge);
                        }
                        else
                        {
                            Triangle &second = mesh.get_triangle(critical_edge_tetra.second);
                            for(int i=0; i<3; i++)
                            {
                                if(!(second.has_vertex(first.TV(i))))
                                {
                                    first.TE(i,critical_edge);
                                    break;
                                }
                            }
                        }
                        //filt_i=filtration[critical_edge[0]-1];
                       // filt_j=filtration[critical_edge[1]-1];
                        filt_i=(filtration[critical_edge[0]-1]>filtration[critical_edge[1]-1])?filtration[critical_edge[0]-1]:filtration[critical_edge[1]-1];
                        filt_j=(filtration[critical_edge[0]-1]>filtration[critical_edge[1]-1])?filtration[critical_edge[1]-1]:filtration[critical_edge[0]-1];
                        Triangle t=mesh.get_triangle(arco->getNode_j()->get_critical_index());
                        for(int i=0;i<3;i++)
                            filt_ex.push_back(filtration[t.TV(i)-1]);
                         index_j=get_max_elevation_vertex(t);
                         
                         
                         
                         }
                    
                    val=abs(mesh.get_vertex(index_i).get_z()-mesh.get_vertex(index_j).get_z());
                    
                    ///////COMMENTED FOR DEBUG
                    if(val <= persistence) /// NEW <= instead of < (uniform execution pattern)
                    {
                        /// this must be enable if we simplify only topologically!! (the same holds in the removal function)
                        Topo_Sempl ts = Topo_Sempl(arco, val, !is_minimum,filt_i,filt_j,filt_ex);
                        q.push(ts);
                    }

                }
            }
            else
            {
                existing_arc->setLabel(2);
            }

        }

        /// set again to begin
        it=extrema->begin();
    }
}

/// push into the priority queue only those arc below the avg
void Forman_Gradient_Simplifier::build_persistence_queue(priority_arcs_queue &q, leaf_ET &local_et, Mesh &mesh, Forman_Gradient &gradient, coord_type persistence)
{
    for(leaf_ET::iterator it_e=local_et.begin(); it_e!=local_et.end(); ++it_e)
    {
        if(gradient.is_edge_critical(it_e->first,it_e->second,mesh))
        {
//            cout<<"saddle --> edge: "<<it_e->first[0]<<" "<<it_e->first[1];
//            cout<<" -- et: "<<it_e->second.first<<" "<<it_e->second.second<<endl;
            /// we get the internal saddle
            iNode* saddle;
            if(it_e->second.second != -1)
                saddle = forman_ig.find_saddle(it_e->second);
            else
            {
                Triangle& top_first = mesh.get_triangle(it_e->second.first);
                short e_pos_first = top_first.edge_index(it_e->first);
                pair<itype,itype> p = make_pair(it_e->second.first, -e_pos_first - 1);
                saddle = forman_ig.find_saddle(p);
            }

            /// and the arcs that link it with minima and maxima
            if(saddle!=NULL)
            {
                coord_type val;

                set<Arc*> &arc_down = saddle->getArcs(false); /// we get the arcs from the saddle to the connected maxima
                for(set<Arc*>::iterator it = arc_down.begin(); it!=arc_down.end(); ++it)
                {
                    if((*it)->getLabel() == 1)
                    {
                        Vertex &v1 = mesh.get_vertex((*it)->getNode_i()->get_critical_index());
                      //  ivect t;
                       // mesh.get_triangle((*it)->getNode_j()->get_critical_index()).convert_to_vec(t);
                        Vertex &v2 = mesh.get_vertex(get_max_elevation_vertex(mesh.get_triangle((*it)->getNode_j()->get_critical_index())));
                        val = abs(v1.get_z() - v2.get_z());

                        if(val <= persistence) /// NEW <= instead of <
                        {
                            Triangle t=mesh.get_triangle((*it)->getNode_j()->get_critical_index());
                          //  cout<<"saddle index: "<<(*it)->getNode_i()->get_critical_index()<<endl;
                         //   cout<<"Triangle index: "<<t.TV(0)<<"; "<<t.TV(1)<<"; "<<t.TV(2)<<endl;
                        int filt0=(filtration[it_e->first[0]-1]>filtration[it_e->first[1]-1])?filtration[it_e->first[0]-1]:filtration[it_e->first[1]-1];
                        int filt1=(filtration[it_e->first[0]-1]>filtration[it_e->first[1]-1])?filtration[it_e->first[1]-1]:filtration[it_e->first[0]-1];
                        ivect filt_ex;
                        

                        for (int i=0;i<3;i++)
                            {filt_ex.push_back(filtration[t.TV(i)-1]);
                         // cout<<"Filt ex "<<i<<" :"<<filtration[t.TV(i)-1]<<endl;
                            }
                        sort(filt_ex.begin(), filt_ex.end(), greater<int>()); 
                        Topo_Sempl ts = Topo_Sempl(*it, val, 1,filt0,filt1,filt_ex);
                          // cout<<"SADDLE is "<< *saddle<<endl;
                          //  cout<<"persistence value is:"<<val<<endl;
                            q.push(ts);
                        }
                    }
                }


                set<Arc*> &arc_up = saddle->getArcs(true); /// we get the arcs from the connected minima to the saddle
                for(set<Arc*>::iterator it = arc_up.begin(); it!=arc_up.end(); ++it)
                {
                    if((*it)->getLabel() == 1)
                    {
                        Vertex &v1 = mesh.get_vertex((*it)->getNode_i()->get_critical_index());
                        Vertex &v2 = mesh.get_vertex((*it)->getNode_j()->get_critical_index());
                        val = abs(v1.get_z() - v2.get_z());
     
                        if(val <= persistence) /// NEW <= instead of <
                        {
                        int filt0=(filtration[it_e->first[0]-1]>filtration[it_e->first[1]-1])?filtration[it_e->first[0]-1]:filtration[it_e->first[1]-1];
                        int filt1=(filtration[it_e->first[0]-1]>filtration[it_e->first[1]-1])?filtration[it_e->first[1]-1]:filtration[it_e->first[0]-1];
                            ivect filt_ex;
                            filt_ex.push_back(filtration[(*it)->getNode_i()->get_critical_index()-1]);

                            Topo_Sempl ts = Topo_Sempl(*it, val, 0,filt0,filt1,filt_ex);
                            // cout<<"SADDLE is "<< *saddle<<endl;
                            // cout<<"persistence value is:"<<val<<endl;
                            q.push(ts);
                        }
                    }
                }

            }
            else
            {
                cout<<"[build_persistence_queue] something wrong here...."<<endl;
                cout<<"missing saddle --> edge: "<<it_e->first[0]<<" "<<it_e->first[1];
                cout<<" -- et: "<<it_e->second.first<<" "<<it_e->second.second<<endl;
                cout<<"et0: "<<mesh.get_triangle(it_e->second.first)<<endl;
                cout<<"et1: "<<mesh.get_triangle(it_e->second.second)<<endl;
                pair<itype,itype> opposite = make_pair(it_e->second.second,it_e->second.first);
                if(forman_ig.find_saddle(opposite)==NULL)
                    cout<<"     the opposite is not contained..."<<endl;
                else
                    cout<<"     the opposite is contained into the MIG"<<endl;
                int a; cin>>a;
            }
        }
    }
}

