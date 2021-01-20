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

void Contraction_Simplifier::find_candidate_edges(Node_V &n, Mesh &mesh, leaf_VT &vts, edge_queue &edges, contraction_parameters &params)
{
    map<ivect, coord_type> lengths;
    ivect e;

    for (RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const &t_id = itPair.first;

        if (mesh.is_triangle_removed(*t_id))
        {
            //  cout<<"triangle removed"<<endl;
            continue;
        }
        Triangle &t = mesh.get_triangle(*t_id);

        for (int i = 0; i < t.vertices_num(); i++)
        {
            t.TE(i, e);

            if (n.indexes_vertex(e[1])) // e (v1,v2) is a candidate edge if at least v2 is in n
            {
                map<ivect, coord_type>::iterator it = lengths.find(e);
                // cout<<e[0]<<" and "<<e[1]<<endl;
                if (it == lengths.end())
                {

                    coord_type length;
                    Vertex &v1 = mesh.get_vertex(e[0]);
                    Vertex &v2 = mesh.get_vertex(e[1]);
                    dvect dif = {v1.get_x() - v2.get_x(), v1.get_y() - v2.get_y(), v1.get_z() - v2.get_z()};
                    //  cout<<dif[0]<<", "<<dif[1]<<", "<<dif[2]<<endl;
                    length = sqrt(dif[0] * dif[0] + dif[1] * dif[1] + dif[2] * dif[2]);

                    //  Edge e((*it)[0],(*it)[1]);
                    lengths[e] = length;
                    //Edge edge_obj(e[0],e[1]);
                    if (length - params.get_maximum_limit() < Zero)
                    {
                        // cout<<"["<<e[0]<<","<<e[1]<<"]  Edge length: "<<length<<endl;
                        edges.push(new Geom_Edge(e, length));
                        //     cout<<"ENQUEUE"<<endl;
                    }
                }
            }
        }
    }
    //  cout<<"======NEXT NODE======"<<endl;
}

void Contraction_Simplifier::find_candidate_edges_parallel(Node_V &n, Mesh &mesh, leaf_VT &vts, edge_queue &edges, contraction_parameters &params, bool is_cross)
{
    map<ivect, coord_type> lengths;
    ivect e;
   /* if (!is_cross)
    {
        for (RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {

            RunIterator const &t_id = itPair.first;

            Triangle &t = mesh.get_triangle(*t_id);
            if (n.partial_indexes_triangle_vertices(t))
            {
                {
                    omp_set_lock(&(t_locks[*t_id - 1]));
                    if (mesh.is_triangle_removed(*t_id))
                    {
                        omp_unset_lock(&(t_locks[*t_id - 1]));
                        continue;
                    }
                    else if (!mesh.is_triangle_removed(*t_id))
                    {
                        for (int i = 0; i < 3; i++)
                        {
                            t.TE(i, e);
                            if (n.completely_indexes_simplex(e)) // e (v1,v2) is a candidate edge if at least v2 is in n
                            {
                                map<ivect, coord_type>::iterator it = lengths.find(e);
                                if (it == lengths.end())
                                {
                                    coord_type length;
                                    Vertex &v1 = mesh.get_vertex(e[0]);
                                    Vertex &v2 = mesh.get_vertex(e[1]);
                                    dvect dif = {v1.get_x() - v2.get_x(), v1.get_y() - v2.get_y(), v1.get_z() - v2.get_z()};
                                    //  cout<<dif[0]<<", "<<dif[1]<<", "<<dif[2]<<endl;
                                    length = sqrt(dif[0] * dif[0] + dif[1] * dif[1] + dif[2] * dif[2]);
                                    lengths[e] = length;

                                    if (length - params.get_maximum_limit() < Zero)
                                    {
                                        // cout<<"["<<e[0]<<","<<e[1]<<"]  Edge length: "<<length<<endl;
                                        edges.push(new Geom_Edge(e, length));
                                        //     cout<<"ENQUEUE"<<endl;
                                    }
                                }
                            }
                        }
                    }
                    omp_unset_lock(&(t_locks[*t_id - 1]));
                }
            }
            else if (n.completely_indexes_triangle_vertices(t))
            { // then there is  no need to check if the edge is contained by the node

                if (mesh.is_triangle_removed(*t_id))
                {
                    continue;
                }

                for (int i = 0; i < 3; i++)
                {
                    t.TE(i, e);

                    map<ivect, coord_type>::iterator it = lengths.find(e);
                    // cout<<e[0]<<" and "<<e[1]<<endl;
                    if (it == lengths.end())
                    {
                        coord_type length;
                        Vertex &v1 = mesh.get_vertex(e[0]);
                        Vertex &v2 = mesh.get_vertex(e[1]);
                        dvect dif = {v1.get_x() - v2.get_x(), v1.get_y() - v2.get_y(), v1.get_z() - v2.get_z()};
                        //  cout<<dif[0]<<", "<<dif[1]<<", "<<dif[2]<<endl;
                        length = sqrt(dif[0] * dif[0] + dif[1] * dif[1] + dif[2] * dif[2]);
                        //  Edge e((*it)[0],(*it)[1]);
                        lengths[e] = length;
                        //Edge edge_obj(e[0],e[1]);
                        if (length - params.get_maximum_limit() < Zero)
                        {
                            // cout<<"["<<e[0]<<","<<e[1]<<"]  Edge length: "<<length<<endl;
                            edges.push(new Geom_Edge(e, length));
                            //     cout<<"ENQUEUE"<<endl;
                        }
                    }
                }
            }
        }
    }
    else
    { */
        for (RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const &t_id = itPair.first;

            if (mesh.is_triangle_removed(*t_id))
            {
                //  cout<<"triangle removed"<<endl;
                continue;
            }
            Triangle &t = mesh.get_triangle(*t_id);

            for (int i = 0; i < t.vertices_num(); i++)
            {
                t.TE(i, e);

                if (n.indexes_vertex(e[1])) // e (v1,v2) is a candidate edge if at least v2 is in n
                {
                    map<ivect, coord_type>::iterator it = lengths.find(e);
                    // cout<<e[0]<<" and "<<e[1]<<endl;
                    if (it == lengths.end())
                    {

                        coord_type length;
                        Vertex &v1 = mesh.get_vertex(e[0]);
                        Vertex &v2 = mesh.get_vertex(e[1]);
                        dvect dif = {v1.get_x() - v2.get_x(), v1.get_y() - v2.get_y(), v1.get_z() - v2.get_z()};
                        //  cout<<dif[0]<<", "<<dif[1]<<", "<<dif[2]<<endl;
                        length = sqrt(dif[0] * dif[0] + dif[1] * dif[1] + dif[2] * dif[2]);

                        //  Edge e((*it)[0],(*it)[1]);
                        lengths[e] = length;
                        //Edge edge_obj(e[0],e[1]);
                        if (length - params.get_maximum_limit() < Zero)
                        {
                            // cout<<"["<<e[0]<<","<<e[1]<<"]  Edge length: "<<length<<endl;
                            edges.push(new Geom_Edge(e, length));
                            //     cout<<"ENQUEUE"<<endl;
                        }
                    }
                }
            }
        }
    //}

    //  cout<<"= =====NEXT NODE======"<<endl;
}

void Contraction_Simplifier::contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1, Node_V &outer_v_block, edge_queue &edges,
                                           Node_V &n, Mesh &mesh, contraction_parameters &params)
{
    //   cout<<"[EDGE CONTRACTION] v1 and v2:"<<e[0]-1<<", "<<e[1]-1<<endl;
    // cout<<"[NOTICE] Contract Edge"<<endl;
    omp_set_lock(&(v_locks[e[0] - 1]));
    omp_set_lock(&(v_locks[e[1] - 1]));
    ivect et_vec;
    et_vec.push_back(et.first);
    if (et.second != -1)
        et_vec.push_back(et.second);

    difference_of_vectors(vt0, et_vec); // vt0 now contains the difference VT0 - ET
    difference_of_vectors(vt1, et_vec); // vt1 now contains the difference VT1 - ET

    // cout<<"VT1 size: "<<vt0.size()<<" VT2 size: "<<vt1.size()<<endl;
    // contract v1 to v0.

    /// prior checking the d-1 faces we update
    /// (1) the corresponding vt0 relation (by adding the triangles in vt1-et to vt0)
    /// (2) then the outer_v_block (if e is a cross edge)
    /// (3) and, finally, we add the new edges crossing b to the edge queue
    if (!params.is_parallel())
        Contraction_Simplifier::update(e, vt0, vt1, n, outer_v_block, edges, mesh, params);
    else
    {

        Contraction_Simplifier::update_parallel(e, vt0, vt1, n, outer_v_block, edges, mesh, params);
    }

    // we remove v2 and the triangles in et
    Contraction_Simplifier::remove_from_mesh(e[1], et, mesh, params);
    // finally we clear the VT(v2)
    vt1.clear();
    //et.clear();
    omp_unset_lock(&(v_locks[e[0] - 1]));
    omp_unset_lock(&(v_locks[e[1] - 1]));
}

void Contraction_Simplifier::get_ET(ivect &e, ET &et, Node_V &n, Mesh &mesh, leaf_VT &vts)
{
    int other_v, local_v_id;
    if (n.indexes_vertex(e[0]))
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
    //  #pragma omp critical
    {
        for (unsigned i = 0; i < vt.size(); i++)
        {

            Triangle &t = mesh.get_triangle(vt[i]);
            if (t.has_vertex(other_v) && (!mesh.is_triangle_removed(vt[i])))
                et_tmp.push_back(vt[i]);
        }
    }

    if (et_tmp.size() == 2)
        et = make_pair(et_tmp[0], et_tmp[1]);
    else if (et_tmp.size() == 0) // Can be deleted, since the default is (-1,-1)
    {
        et = make_pair(-1, -1);
    }
    else
    {
        et = make_pair(et_tmp[0], -1);
    }
}

void Contraction_Simplifier::clean_coboundary(VT &cob, Mesh &mesh)
{

    for (ivect_iter it = cob.begin(); it != cob.end();)
    {
        int tid = *it;
        //     omp_set_lock(&(t_locks[tid - 1]));  // It is possible that the triangle is not in either n or outer_v_block
        if (mesh.is_triangle_removed(*it))
            cob.erase(it);
        else
            ++it;
        //   omp_unset_lock(&(t_locks[tid - 1]));
    }
}

void Contraction_Simplifier::update(const ivect &e, VT &vt, VT &difference, Node_V &n, Node_V &v_block,
                                    edge_queue &edges, Mesh &mesh, contraction_parameters &params)
{
    set<ivect> e_set; /// we insert the new edges first in this set to avoid duplicate insertions in the queue

    if (params.is_QEM())
    {
        initialQuadric[e[0]] += initialQuadric[e[1]]; // Update the QEM of e[0]
        set<ivect> v1_related_set;
        for (int i = 0; i < vt.size(); i++)
        { //Store VE(v1)
            ivect new_e(2, 0);
            Triangle &t = mesh.get_triangle(vt[i]);
            int v1_pos = t.vertex_index(e[0]);
            for (int i = 0; i < t.vertices_num(); i++)
            {
                if (i != v1_pos)
                {
                    t.TE(i, new_e);
                    //     if(n.indexes_vertex(new_e[1])) /// we process an edge only if it has all the extrema already processed
                    v1_related_set.insert(new_e);
                }
            }
        }
        for (auto it = v1_related_set.begin(); it != v1_related_set.end(); ++it)
        {
            //Calculate updated edge costs of VE(v1)
            double value;
            ivect e(2, -1);
            Vertex &v1 = mesh.get_vertex((*it)[0]);
            Vertex &v2 = mesh.get_vertex((*it)[1]);
            int new_vertex_pos = -1;
            double error = compute_error((*it)[0], (*it)[1], mesh, new_vertex_pos);
            //    cout<<"[DEBUG] calculated error: "<<error<<endl;
            assert(new_vertex_pos != -1);

            if (updated_edges.find(*it) != updated_edges.end())
                updated_edges[*it] = error;
            else
            {
                pair<ivect, double> updated_edge(*it, error);
                updated_edges.insert(updated_edge);
            }

            e = {(*it)[0], (*it)[1]};
            if ((error - params.get_maximum_limit() < Zero) && n.indexes_vertex(e[1]))
            {
                if (new_vertex_pos == 1)
                {
                    e = {(*it)[1], (*it)[0]};
                }

                //  cout<<"["<<e[0]-1<<","<<e[1]-1<<"]  Error will be introduced: "<<error<<endl;
                edges.push(new Geom_Edge(e, error));
            }
        }
    }
    //  #pragma omp critical
    {
        for (ivect_iter it = difference.begin(); it != difference.end(); ++it)
        {
            int pos = -1;
            //   bool not_valid=false;
            Triangle &t = mesh.get_triangle(*it);
            if (!n.completely_indexes_triangle_vertices(t))
            {
                //#pragma omp critical
                {
                    //   omp_set_lock(&(t_locks[*it - 1]));
                    if (!mesh.is_triangle_removed(*it))
                    {
                        /// then we update the triangle changing e[1] with e[0]
                        pos = t.vertex_index(e[1]);
                        t.setTV_keep_border(pos, e[0]);
                        //       omp_unset_lock(&(t_locks[*it - 1]));
                    }
                    else
                    {
                        //        omp_unset_lock(&(t_locks[*it - 1]));
                        continue;
                    }
                }
                // if(not_valid)
                //     continue;
            }
            else
            {
                /// then we update the triangle changing e[1] with e[0]
                pos = t.vertex_index(e[1]);
                t.setTV_keep_border(pos, e[0]);
            }
            /// before updating the triangle, we check
            /// if the leaf block indexing e[0] does not contain the current triangle we have to add it
            /// NOTA: there is one possible case.. as leaf block n already indexes the triangle in e[1]
            /// NOTA2: we have just to check that n and v_block are different (as if they are equal the edge is internal in n)
            // Note that in the triangle *it, the e[1] has been replaced by e[0] already

            if (n.get_v_start() != v_block.get_v_start())
            {
                if (e[0] < e[1] && !v_block.index_triangle(*it, mesh, e[0]))
                {
                    v_block.add_triangle(*it);
                }
                else if (e[0] > e[1] && !n.index_triangle(*it, mesh, e[0]))
                {
                    n.add_triangle(*it);
                }
            }

            dvect diff(4, 0.0);
            if (params.is_QEM())
            {
                ivect new_e(2, 0); // new_e.assign(2,0);
                for (int i = 0; i < 3; i++)
                {
                    if (i != pos)
                    {
                        mesh.get_triangle(*it).TE(i, new_e); //t.TE(new_e,pos,i); need to check
                        e_set.insert(new_e);
                    }
                }
            }
            else
            {
                /// we have to add the new edges in the queue
                ivect new_e;
                new_e.assign(2, 0);
                if (!n.completely_indexes_triangle_vertices(mesh.get_triangle(*it)))
                {
                    //#pragma omp critical
                    // omp_set_lock(&(t_locks[*it - 1]));
                    {
                        for (int i = 0; i < 3; i++)
                        {
                            if (i != pos)
                            {
                                if (!mesh.get_triangle(*it).is_removed())
                                {
                                    mesh.get_triangle(*it).TE(i, new_e); //t.TE(new_e,pos,i); need to check
                                    e_set.insert(new_e);
                                }
                            }
                        }
                    }
                    //   omp_unset_lock(&(t_locks[*it - 1]));
                }
                else
                {
                    for (int i = 0; i < 3; i++)
                    {
                        if (i != pos)
                        {
                            mesh.get_triangle(*it).TE(i, new_e); //t.TE(new_e,pos,i); need to check
                            e_set.insert(new_e);
                        }
                    }
                }
            }
        }
    }

    /// we push the new "unique" edges in the queue
    for (auto it = e_set.begin(); it != e_set.end(); ++it)
    {
        //Calculate length
        double value;
        ivect e(2, -1);

        if (params.is_QEM())
        {
            int new_vertex_pos = -1;
            double error = compute_error((*it)[0], (*it)[1], mesh, new_vertex_pos);
            //    cout<<"[DEBUG] calculated error: "<<error<<endl;
            assert(new_vertex_pos != -1);
            if (updated_edges.find(*it) != updated_edges.end()) // updated_edges keep all the edges that have new cost values
                updated_edges[*it] = error;
            else
            {
                pair<ivect, double> updated_edge(*it, error);
                updated_edges.insert(updated_edge);
            }
            e = {(*it)[0], (*it)[1]};
            if ((error - params.get_maximum_limit() < Zero) && n.indexes_vertex(e[1]))
            {
                if (new_vertex_pos == 1)
                {
                    e = {(*it)[1], (*it)[0]};
                }

                //  cout<<"["<<e[0]-1<<","<<e[1]-1<<"]  Error will be introduced: "<<error<<endl;

                edges.push(new Geom_Edge(e, error));
            }
        }
        else
        {
            if (!params.is_parallel())
            {
                Vertex &v1 = mesh.get_vertex((*it)[0]);
                Vertex &v2 = mesh.get_vertex((*it)[1]);
                dvect dif = {v1.get_x() - v2.get_x(), v1.get_y() - v2.get_y(), v1.get_z() - v2.get_z()};
                value = sqrt(dif[0] * dif[0] + dif[1] * dif[1] + dif[2] * dif[2]);
                e = {(*it)[0], (*it)[1]};

                if ((value - params.get_maximum_limit() < Zero) && n.indexes_vertex(e[1]))
                {
                    // Geom_Edge new_edge(e,length);
                    cout << "[" << e[0] << "," << e[1] << "]  Edge length: " << value << endl;
                    edges.push(new Geom_Edge(e, value));
                }
            }
            else
            {
                if (n.completely_indexes_simplex((*it)))
                {
                    Vertex &v1 = mesh.get_vertex((*it)[0]);
                    Vertex &v2 = mesh.get_vertex((*it)[1]);
                    dvect dif = {v1.get_x() - v2.get_x(), v1.get_y() - v2.get_y(), v1.get_z() - v2.get_z()};
                    value = sqrt(dif[0] * dif[0] + dif[1] * dif[1] + dif[2] * dif[2]);
                    e = {(*it)[0], (*it)[1]};
                    if ((value - params.get_maximum_limit() < Zero))
                    {
                        //     cout<<"["<<e[0]<<","<<e[1]<<"]  Edge length: "<<value<<endl;
                        edges.push(new Geom_Edge(e, value));
                    }
                }
            }
        }
    }

    /// finally we update the VT relation of e[0]
    unify_vectors(vt, difference);
}

void Contraction_Simplifier::remove_from_mesh(int to_delete_v, ET &et, Mesh &mesh, contraction_parameters &params)
{
    if (et.first != -1)
    {
        //omp_set_lock(&(t_locks[et.first - 1]));
        mesh.remove_triangle(et.first);
      //  omp_unset_lock(&(t_locks[et.first - 1]));
        params.increment_counter();
    }
    if (et.second != -1)
    {
       // omp_set_lock(&(t_locks[et.second - 1]));
        mesh.remove_triangle(et.second);
      //  omp_unset_lock(&(t_locks[et.second - 1]));
        params.increment_counter();
    }

    mesh.remove_vertex(to_delete_v);

    params.increment_contracted_edges_counter();
}

bool Contraction_Simplifier::link_condition(int v0, int v1, VT &vt0, VT &vt1, ET &et, Mesh &mesh)
{
    //iset link_ab;
    //  iset link_e;
    //iset vts;

    //Update: Considering that the edge to be contracted should not be boundary edge
    //We can simplify the link condition check while checking if e is boundary edge
    if (et.first == -1 || et.second == -1)
        return false;
    int counter = 0;
    //#pragma omp critical
    {
        iset vv0, vv1;
        //vts.insert(vts.end(),vt0.begin(),vt0.end());
        VT vts = vt0;
        unify_vectors(vts, vt1);
        for (int i = 0; i < vts.size(); i++)
        {
            omp_set_lock(&(t_locks[vts[i] - 1]));
        }

        for (int i = 0; i < vt0.size(); i++)
        {
            if (!mesh.is_triangle_removed(vt0[i]))
            {
                Triangle t = mesh.get_triangle(vt0[i]);
                int v0_id = t.vertex_index(v0);
                int v01_id = t.TV((v0_id + 1) % 3);
                int v02_id = t.TV((v0_id + 2) % 3);
                vv0.insert(v01_id);
                vv0.insert(v02_id);
            }
        }

        for (int i = 0; i < vt1.size(); i++)
        {
            if (!mesh.is_triangle_removed(vt1[i]))
            {
                Triangle t = mesh.get_triangle(vt1[i]);
                int v1_id = t.vertex_index(v1);

                vv1.insert(t.TV((v1_id + 1) % 3));
                vv1.insert(t.TV((v1_id + 2) % 3));
            }
        }

        //  cout<<v1<<"'s VV size: "<<vv1.size()<<endl;
        for (iset_iter it = vv1.begin(); it != vv1.end(); it++)
        {
            if (vv0.find(*it) != vv0.end())
            {
                //   link_ab.insert(*it);
                counter++;
            }
        }

        for (int i = 0; i < vts.size(); i++)
        {
            omp_unset_lock(&(t_locks[vts[i] - 1]));
        }

        // for (int i=0;i<vt1.size();i++){
        //     omp_unset_lock(&(t_locks[vt1[i]]));
        // }
        // if(et.first!=-1){
        //     if(!mesh.is_triangle_removed(et.first)){
        //     Triangle t1= mesh.get_triangle(et.first);
        //     for(int i=0;i<3;i++){
        //         if(t1.TV(i)!=v0&&t1.TV(i)!=v1)
        //          link_e.insert(t1.TV(i));
        //     }
        //     }
        // }
        // if(et.second!=-1){
        //      if(!mesh.is_triangle_removed(et.second)){
        //     Triangle t2= mesh.get_triangle(et.second);
        //     for(int i=0;i<3;i++){
        //         if(t2.TV(i)!=v0&&t2.TV(i)!=v1)
        //          link_e.insert(t2.TV(i));
        //     }
        // }
        // }
    }
    return counter <= 2;
}

bool Contraction_Simplifier::link_condition(int v0, int v1, VT &vt0, VT &vt1, ET &et, Node_V &n, VV &vv_locks, Mesh &mesh)
{

    //Update: Considering that the edge to be contracted should not be boundary edge
    //We can simplify the link condition check while checking if e is boundary edge
    if (et.first == -1 || et.second == -1)
        return false;
    int counter = 0;
    //#pragma omp critical
    {
        iset vv0, vv1;
        //vts.insert(vts.end(),vt0.begin(),vt0.end());
        VT vts = vt0;
        unify_vectors(vts, vt1);
        // for (int i = 0; i < vts.size(); i++)
        // {
        //     omp_set_lock(&(t_locks[vts[i] - 1]));
        // }

        for (int i = 0; i < vt0.size(); i++)
        {
            if (!mesh.is_triangle_removed(vt0[i]))
            {
                Triangle t = mesh.get_triangle(vt0[i]);
                int v0_id = t.vertex_index(v0);
                int v01_id = t.TV((v0_id + 1) % 3);
                int v02_id = t.TV((v0_id + 2) % 3);
                if (!n.indexes_vertex(v01_id))
                    vv_locks.insert(v01_id);

                if (!n.indexes_vertex(v02_id))
                    vv_locks.insert(v02_id);

                vv0.insert(v01_id);
                vv0.insert(v02_id);
            }
        }

        for (int i = 0; i < vt1.size(); i++)
        {
            if (!mesh.is_triangle_removed(vt1[i]))
            {
                Triangle t = mesh.get_triangle(vt1[i]);
                int v1_id = t.vertex_index(v1);
                int v11_id = t.TV((v1_id + 1) % 3);
                int v12_id = t.TV((v1_id + 2) % 3);
                if (!n.indexes_vertex(v11_id))
                    vv_locks.insert(v11_id);
                if (!n.indexes_vertex(v12_id))
                    vv_locks.insert(v12_id);
                vv1.insert(v11_id);
                vv1.insert(v12_id);
            }
        }

        //  cout<<v1<<"'s VV size: "<<vv1.size()<<endl;
        for (iset_iter it = vv1.begin(); it != vv1.end(); it++)
        {
            if (vv0.find(*it) != vv0.end())
            {
                //   link_ab.insert(*it);
                counter++;
            }
        }

        // for (int i = 0; i < vts.size(); i++)
        // {
        //     omp_unset_lock(&(t_locks[vts[i] - 1]));
        // }

        for (iset_iter it = vv_locks.begin(); it != vv_locks.end(); it++)
        {
            omp_set_lock(&(v_locks[*it - 1]));
        }

        // for (int i=0;i<vt1.size();i++){
        //     omp_unset_lock(&(t_locks[vt1[i]]));
        // }
        // if(et.first!=-1){
        //     if(!mesh.is_triangle_removed(et.first)){
        //     Triangle t1= mesh.get_triangle(et.first);
        //     for(int i=0;i<3;i++){
        //         if(t1.TV(i)!=v0&&t1.TV(i)!=v1)
        //          link_e.insert(t1.TV(i));
        //     }
        //     }
        // }
        // if(et.second!=-1){
        //      if(!mesh.is_triangle_removed(et.second)){
        //     Triangle t2= mesh.get_triangle(et.second);
        //     for(int i=0;i<3;i++){
        //         if(t2.TV(i)!=v0&&t2.TV(i)!=v1)
        //          link_e.insert(t2.TV(i));
        //     }
        // }
        // }
    }
    return counter <= 2;
}

bool Contraction_Simplifier::link_condition(int v0, int v1, VT &vt0, VT &vt1, ET &et, Node_V &n, Node_V &v_block, VV &vv_locks, Mesh &mesh)
{

    //Update: Considering that the edge to be contracted should not be boundary edge
    //We can simplify the link condition check while checking if e is boundary edge
    if (et.first == -1 || et.second == -1)
        return false;
    int counter = 0;
    //#pragma omp critical
    //  cout<<"enter link condition function"<<endl;
    {
        iset vv0, vv1;
        //vts.insert(vts.end(),vt0.begin(),vt0.end());
        VT vts = vt0;

        unify_vectors(vts, vt1);
        // for (int i = 0; i < vts.size(); i++)
        // {
        //     omp_set_lock(&(t_locks[vts[i] - 1]));
        // }

        for (int i = 0; i < vt0.size(); i++)
        {
            if (!mesh.is_triangle_removed(vt0[i]))
            {
                Triangle t = mesh.get_triangle(vt0[i]);
                int v0_id = t.vertex_index(v0);
                int v01_id = t.TV((v0_id + 1) % 3);
                int v02_id = t.TV((v0_id + 2) % 3);

                if (!n.indexes_vertex(v01_id) && !v_block.indexes_vertex(v01_id))
                    vv_locks.insert(v01_id);

                if (!n.indexes_vertex(v02_id) && !v_block.indexes_vertex(v02_id))
                    vv_locks.insert(v02_id);

                vv0.insert(v01_id);
                vv0.insert(v02_id);
            }
        }

        for (int i = 0; i < vt1.size(); i++)
        {
            if (!mesh.is_triangle_removed(vt1[i]))
            {
                Triangle t = mesh.get_triangle(vt1[i]);
                int v1_id = t.vertex_index(v1);
                int v11_id = t.TV((v1_id + 1) % 3);
                int v12_id = t.TV((v1_id + 2) % 3);
                if (!n.indexes_vertex(v11_id) && !v_block.indexes_vertex(v11_id))
                    vv_locks.insert(v11_id);
                if (!n.indexes_vertex(v12_id) && !v_block.indexes_vertex(v12_id))
                    vv_locks.insert(v12_id);
                vv1.insert(v11_id);
                vv1.insert(v12_id);
            }
        }

        //  cout<<v1<<"'s VV size: "<<vv1.size()<<endl;
        for (iset_iter it = vv1.begin(); it != vv1.end(); it++)
        {
            if (vv0.find(*it) != vv0.end())
            {
                //   link_ab.insert(*it);
                counter++;
            }
        }

        for (iset_iter it = vv_locks.begin(); it != vv_locks.end(); it++)
        {
            omp_set_lock(&(v_locks[*it - 1]));
        }
    }
    return counter <= 2;
}

void Contraction_Simplifier::update_parallel(const ivect &e, VT &vt, VT &difference, Node_V &n, Node_V &v_block, edge_queue &edges,
                                             Mesh &mesh, contraction_parameters &params)
{
    set<ivect> e_set; /// we insert the new edges first in this set to avoid duplicate insertions in the queue
    ivect t_locks_id;
    /// before updating the triangle, we check
    /// if the leaf block indexing e[0] does not contain the current triangle we have to add it
    /// NOTA: there is one possible case.. as leaf block n already indexes the triangle in e[1]
    /// NOTA2: we have just to check that n and v_block are different (as if they are equal the edge is internal in n)

    for (ivect_iter it = difference.begin(); it != difference.end(); ++it)
    {
        int pos = -1;
        //   bool not_valid=false;
        dvect diff(4, 0.0);
        /// we have to add the new edges in the queue
        ivect new_e;
        new_e.assign(2, 0);

        Triangle &t = mesh.get_triangle(*it);
        if (n.get_v_start() != v_block.get_v_start())
        {
            if (e[0] < e[1] && !v_block.index_triangle(*it, mesh, e[1]))
            {
                v_block.add_triangle(*it);
            }
            else if (e[0] > e[1] && !n.index_triangle(*it, mesh, e[1]))
            {
                n.add_triangle(*it);
            }
        }

        if (!mesh.is_triangle_removed(*it))
        {
            /// then we update the triangle changing e[1] with e[0]
            pos = t.vertex_index(e[1]);
            t.setTV_keep_border(pos, e[0]);
            for (int i = 0; i < 3; i++)
            {
                if (i != pos)
                {
                    t.TE(i, new_e); //t.TE(new_e,pos,i); need to check
                    e_set.insert(new_e);
                }
            }
        }
        else
        {
            //    omp_unset_lock(&(t_locks[*it - 1]));
            continue;
        }
    }

    /// we push the new "unique" edges in the queue
    for (auto it = e_set.begin(); it != e_set.end(); ++it)
    {
        //Calculate length
        double value;
        ivect e(2, -1);
       // if (n.indexes_vertex(e[1]))
        //{
            Vertex &v1 = mesh.get_vertex((*it)[0]);
            Vertex &v2 = mesh.get_vertex((*it)[1]);
            dvect dif = {v1.get_x() - v2.get_x(), v1.get_y() - v2.get_y(), v1.get_z() - v2.get_z()};
            value = sqrt(dif[0] * dif[0] + dif[1] * dif[1] + dif[2] * dif[2]);
            e = {(*it)[0], (*it)[1]};
            if ((value - params.get_maximum_limit() < Zero)&&n.indexes_vertex(e[1]))
            {
                //     cout<<"["<<e[0]<<","<<e[1]<<"]  Edge length: "<<value<<endl;
                edges.push(new Geom_Edge(e, value));
            }
       // }
    }
    // for(ivect_iter it=t_locks_id.begin();it!=t_locks_id.end();it++){
    //     omp_unset_lock(&(t_locks[*it]));
    // }
    /// finally we update the VT relation of e[0]
    unify_vectors(vt, difference);
}

void Contraction_Simplifier::simplify_parallel(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli)
{
   // cerr << "[NOTICED] Cache size: " << cli.cache_size << endl;
   // LRU_Cache<int, leaf_VT> cache(cli.cache_size); // the key is v_start while the value are the VT relations
    contraction_parameters params;
    params.set_maximum_limit(cli.maximum_limit);

    params.queue_criterion_length();
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
//     for (int i = 0; i < t_num; i++)
//         omp_init_lock(&(t_locks[i]));
//     cout << "Initialize t_locks" << endl;
#pragma omp parallel for
    for (int i = 0; i < v_num; i++)
        omp_init_lock(&(v_locks[i]));
    cout << "Initialize v_locks" << endl;

#pragma omp parallel for
    for (int i = 0; i < l_num; i++)
        omp_init_lock(&(l_locks[i]));
    cout << "Initialize l_locks" << endl;

    while (1)
    {
        simplification_round = params.get_contracted_edges_num(); //checked edges
        /// HERE YOU NEED TO DEFINE A PROCEDURE FOR SIMPLIFY THE TIN BY USING THE SPATIAL INDEX
        this->simplify_compute_parallel(mesh,  tree.get_subdivision(), params, tree);

        cout << "Num of edges enqueued:" << params.get_sum_edge_queue_sizes() << endl;
        // PARTIAL SIMPLIFICATION STATS
        cerr << "=== end-of-round " << round << ") --> contracted edges: ";
        cerr << params.get_contracted_edges_num() - simplification_round << endl;
        round++;

        if (cli.debug_mode)
        {
            time.stop();
            time.print_elapsed_time("   [TIME] executing a simplification round: ");
            time.start();
            //  cerr << "   [RAM] peak for executing a simplification round: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
        }
        if (simplification_round == params.get_contracted_edges_num())
            break;


    }

// #pragma omp parallel for
//     for (int i = 0; i < t_num; i++)
//         omp_destroy_lock(&(t_locks[i]));

#pragma omp parallel for
    for (int i = 0; i < v_num; i++)
        omp_destroy_lock(&(v_locks[i]));

#pragma omp parallel for
    for (int i = 0; i < l_num; i++)
        omp_destroy_lock(&(l_locks[i]));

    time.stop();

    ///// Clear all the auxiliary data structures.
    //v_locks.clear();
    vector<omp_lock_t>().swap(v_locks);
    vector<omp_lock_t>().swap(l_locks);
    vector<int>().swap(v_in_leaf);
    lists_leafs().swap(conflict_leafs);
   // l_locks.clear();


    if (!cli.debug_mode)
        time.print_elapsed_time("[TIME] Edge contraction simplification: ");
    cerr << "[MEMORY] peak for Simplification: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
    //else
    //params.print_simplification_partial_timings();
    //params.print_simplification_counters();
    /// finally we have to update/compress the mesh and the tree
    Contraction_Simplifier::update_mesh_and_tree(tree, mesh, params);
    cerr << "[MEMORY] peak for mesh and tree updating: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

}

void Contraction_Simplifier::simplify(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli)
{

    //cerr<<"==Homology preserving simplification - weak-link condition=="<<endl;

    cerr << "[NOTICED] Cache size: " << cli.cache_size << endl;
    LRU_Cache<int, leaf_VT> cache(cli.cache_size); // the key is v_start while the value are the VT relations
    contraction_parameters params;
    params.set_maximum_limit(cli.maximum_limit);
    if (cli.QEM_based)
        params.queue_criterion_QEM();
    else
    {
        params.queue_criterion_length();
    }

    Timer time;
    int simplification_round;
    int round = 1;
    if (params.is_QEM())
    {

        trianglePlane = vector<dvect>(mesh.get_triangles_num(), dvect(4, 0));
        initialQuadric = vector<Matrix>(mesh.get_vertices_num() + 1, Matrix(0.0));
        cout << "=========Calculate triangle plane========" << endl;
        compute_triangle_plane(mesh, trianglePlane);
        cout << "=========Calculate initial QEM========" << endl;
        compute_initial_QEM(mesh, trianglePlane);
    }
    time.start();
    while (1)
    {
        simplification_round = params.get_contracted_edges_num(); //checked edges

        /// HERE YOU NEED TO DEFINE A PROCEDURE FOR SIMPLIFY THE TIN BY USING THE SPATIAL INDEX
        this->simplify_compute(tree.get_root(), mesh, cache, tree.get_subdivision(), params, tree);

        cout << "Num of edges enqueued:" << params.get_sum_edge_queue_sizes() << endl;
        // PARTIAL SIMPLIFICATION STATS
        cerr << "=== end-of-round " << round << ") --> contracted edges: ";
        cerr << params.get_contracted_edges_num() - simplification_round << endl;
        round++;

        if (cli.debug_mode)
        {
            time.stop();
            time.print_elapsed_time("   [TIME] executing a simplification round: ");
            time.start();
            //  cerr << "   [RAM] peak for executing a simplification round: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
        }

        if (simplification_round == params.get_contracted_edges_num())
            break;

        cache.reset();
        //    if(round==2)
        //      break;
    }
    time.stop();
    if (!cli.debug_mode)
        time.print_elapsed_time("[TIME] Edge contraction simplification: ");
    //else
    //params.print_simplification_partial_timings();
    //params.print_simplification_counters();

    //  cerr << "[RAM peak] for contracting a simplicial complex: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

    /// finally we have to update/compress the mesh and the tree
    Contraction_Simplifier::update_mesh_and_tree(tree, mesh, params);
}

void Contraction_Simplifier::simplify_compute(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, Spatial_Subdivision &division, contraction_parameters &params, PRT_Tree &tree)
{
    if (n.is_leaf())
    {
        if (params.is_QEM() == true)
            simplify_leaf_QEM(n, mesh, cache, params, tree);
        else
            simplify_leaf(n, mesh, cache, params, tree);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if (n.get_son(i) != NULL)
            {
                simplify_compute(*n.get_son(i), mesh, cache, division, params, tree);
            }
        }
    }
}

void Contraction_Simplifier::simplify_compute_parallel(Mesh &mesh,  Spatial_Subdivision &division, contraction_parameters &params, PRT_Tree &tree)
{

    //First part for internal edges
    // #pragma omp parallel for // schedule(dynamic,1)
    //     for (unsigned i = 0; i < tree.get_leaves_number(); i++)
    //     {

    //         Node_V *leaf = tree.get_leaf(i);

    //         simplify_leaf(*leaf, mesh, cache, params, tree);

    //     }
    //Second part for cross edges

    //boost::dynamic_bitset<> conflict_nodes(tree.get_leaves_number());
    ivect nodes_status(tree.get_leaves_number(), 0);
    bool processed = false;
    do
    {
        processed = false;
#pragma omp parallel for // schedule(dynamic,1)
        for (unsigned i = 0; i < tree.get_leaves_number(); i++)
        {
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
                omp_unset_lock(&(l_locks[i]));
                iset conflicts = conflict_leafs[i];
                // Check if the conflict nodes were set to 1 already
                bool shared_conflicts = false;
                for (iset_iter it = conflicts.begin(); it != conflicts.end(); it++)
                {
                    //  cout<<"set leaf node:"<<*it<<" on thread "<<omp_get_thread_num()<<endl;
                    omp_set_lock(&(l_locks[*it])); // leafs should be locked when being checked and updated.
                    int status = 0;
                    //   #pragma omp atomic read
                    status = nodes_status[*it];
                    if (status == 1) //it is conflicted with another node being processed
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
                        shared_conflicts = true;
                        //cout<<"neighbor conflict"<<endl;

                        break;
                    }
                    else if (status == 0)
                    {
                        nodes_status[*it] = 1;
                    }
                    omp_unset_lock(&(l_locks[*it]));
                }
                if (shared_conflicts == true)
                {
                    // cout<<"unset leaf node:"<<i<<" on thread "<<omp_get_thread_num()<<endl;
                    //  omp_unset_lock(&(l_locks[i]));
                    continue;
                }
                // // Set the status of conflict nodes to 1
                // for(iset_iter it=conflicts.begin();it!=conflicts.end();it++){
                //     omp_set_lock(&(l_locks[*it]));
                //     if(nodes_status[*it]!=-1)
                //         nodes_status[*it]=1;
                //     // cout<<"unset leaf node:"<<*it<<" on thread "<<omp_get_thread_num()<<endl;
                //      omp_unset_lock(&(l_locks[*it]));
                // }
                Node_V *leaf = tree.get_leaf(i);
                //  cout<<"Node "<<i<<" will be processed."<<endl;
                omp_set_lock(&(l_locks[i]));
                // set nodes_status[i]=2 when node is being processed
                nodes_status[i] = 2;
                omp_unset_lock(&(l_locks[i]));

                processed = true;
               // cout << "Start simplification" << endl;
                simplify_leaf_cross(*leaf, i, mesh, params, tree);
               // cout << "Finish simplification" << endl;
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
                        //    #pragma omp atomic write
                        nodes_status[*it] = 0;
                    }
                    //cout<<"unset leaf node:"<<*it<<" on thread "<<omp_get_thread_num()<<endl;
                    omp_unset_lock(&(l_locks[*it]));
                }
            }
        }
       // cout << "finished one for loop" << endl;

    } while (processed == true);
}

void Contraction_Simplifier::simplify_leaf_QEM(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, PRT_Tree &tree)
{

    if (!n.indexes_vertices())
        return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();
    itype v_range = v_end - v_start;

    //cout<<"Simplification in leaf."<<endl;
    // leaf_VT local_vts(v_range,VT());
    // n.get_VT(local_vts,mesh);
    leaf_VT &local_vts = get_VTS(n, mesh, cache, tree, params);
    // Create a priority queue of candidate edges
    edge_queue edges;
    find_candidate_edges_QEM(n, mesh, local_vts, edges, params);
    int edge_num = edges.size();
    int edges_contracted_leaf = 0;
  //  cout << "Edge number:" << edges.size() << endl;
    params.add_edge_queue_size(edges.size());
    while (!edges.empty())
    {
        Geom_Edge *current = edges.top();
        ivect e = current->edge;
        //    cout<<"Start contraction."<<endl;
        //  cout<<"Edge error:"<<current->val<<endl;

        edges.pop();

        if (mesh.is_vertex_removed(e[0]) || mesh.is_vertex_removed(e[1]))
        {
            //       cout<<"skip current edge."<<endl;
            //   cout<<"[DEBUG] edge not complete: "<<e[0]<<", "<<e[1]<<endl;
            //   cout<<"Vertex removed"<<endl;
            delete current;
            // if(edges_contracted_leaf>edge_num*0.2)
            // cout<<"Num of deleted in a leaf"<<edges_contracted_leaf<<"; 20% of the queue num:"<<edge_num*0.2<<endl;
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
                // cout<<"skip current edge."<<endl;
                // cout<<"[DEBUG] edge: "<<sorted_e[0]<<", "<<sorted_e[1]<<"; updated error: "<<it->second<<"old error: "<<current->val<<endl;
                delete current;
                continue;
            }
        }

        ET et(-1, -1);
        VT *vt0 = NULL, *vt1 = NULL;
        Node_V *outer_v_block = NULL;

        get_edge_relations(e, et, vt0, vt1, outer_v_block, n, mesh, local_vts, cache, params, tree);
        if (link_condition(e[0], e[1], *vt0, *vt1, et, mesh))
        {
            contract_edge(e, et, *vt0, *vt1, *outer_v_block, edges, n, mesh, params);
            edges_contracted_leaf++;
            // break;
        }
      //  cout << "Number of edges remaining:" << edges.size() << endl;
    }

    //update the cache if the vts have been stored already.
    // if(cache.find(v_start) != cache.end()){
    //     cache.update(v_start,local_vts);
    // }
}

void Contraction_Simplifier::simplify_leaf(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, PRT_Tree &tree)
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
    ////// THIS method should be used when cache is used in parallel computing
    leaf_VT &local_vts = get_VTS(n, mesh, cache, tree, params);
    //// Here we use the simple get_VT to check if that is the cause of slow computing.

    // leaf_VT local_vts(v_range,VT());
    // n.get_VT(local_vts,mesh);

    // Create a priority queue of candidate edges
    edge_queue edges;
    if (params.is_parallel())
    {
        find_candidate_edges_parallel(n, mesh, local_vts, edges, params, false);
    }
    else
        find_candidate_edges(n, mesh, local_vts, edges, params);
    int edge_num = edges.size();
    int edges_contracted_leaf = 0;
    // cout<<"Edge number:"<<edges.size()<<endl;
    // cout<<"Number of threads used: "<<omp_get_num_threads()<<endl;
    // cout<<"Current thread id: "<<omp_get_thread_num()<<endl;

    params.add_edge_queue_size(edges.size());
    while (!edges.empty())
    {
        Geom_Edge *current = edges.top();
        ivect e = current->edge;
        //    cout<<"Start contraction."<<endl;
        //  cout<<"Edge Length:"<<current->val<<endl;

        edges.pop();
        // NOTE: Should be atomic when we have cross edges.
        if (mesh.is_vertex_removed(e[0]) || mesh.is_vertex_removed(e[1]))
        {

            //   cout<<"Vertex removed"<<endl;
            // cout<<"skip current edge"<<endl;
            delete current;
            // if(edges_contracted_leaf>edge_num*0.2)
            // cout<<"Num of deleted in a leaf"<<edges_contracted_leaf<<"; 20% of the queue num:"<<edge_num*0.2<<endl;
            continue;
        }

        ET et(-1, -1);
        VT *vt0 = NULL, *vt1 = NULL;
        Node_V *outer_v_block = NULL;

        get_edge_relations(e, et, vt0, vt1, outer_v_block, n, mesh, local_vts, cache, params, tree);

        if (params.is_parallel())
        {
            VV vv_locks;
            if (link_condition(e[0], e[1], *vt0, *vt1, et, n, vv_locks, mesh))
            {
                contract_edge(e, et, *vt0, *vt1, *outer_v_block, edges, n, mesh, params);
                edges_contracted_leaf++;
                // break;
            }
            for (iset_iter it = vv_locks.begin(); it != vv_locks.end(); it++)
            {
                omp_unset_lock(&(v_locks[*it - 1]));
            }
        }
        else
        {
            if (link_condition(e[0], e[1], *vt0, *vt1, et, mesh))
            {
                contract_edge(e, et, *vt0, *vt1, *outer_v_block, edges, n, mesh, params);
                edges_contracted_leaf++;
                // break;
            }
        }
    }

    // leaf_VV vvs;
    // n.get_VV(vvs,mesh);
}

void Contraction_Simplifier::simplify_leaf_cross(Node_V &n, int n_id, Mesh &mesh, contraction_parameters &params, PRT_Tree &tree)
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
    if (params.is_parallel())
    {
        find_candidate_edges_parallel(n, mesh, local_vts, edges, params, true);
    }
    else
        find_candidate_edges(n, mesh, local_vts, edges, params);

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
        get_edge_relations(e, et, vt0, vt1,v1_is_border,v2_is_border, outer_v_block, n, mesh, local_vts, is_v_border, local_cache, params, tree);

        // if(params.is_parallel()){
        VV vv_locks;
        if (link_condition(e[0], e[1], *vt0, *vt1, et, n, *outer_v_block, vv_locks, mesh)&&valid_boundary_condition(e[0],e[1],*vt0,*vt1,et,v1_is_border,v2_is_border,mesh)&&not_fold_over(e[0], e[1], *vt0, *vt1, et, mesh))
        {
            contract_edge(e, et, *vt0, *vt1, *outer_v_block, edges, n, mesh, params);
            edges_contracted_leaf++;
            // break;

            // A new step for cross edge case
            // Check possible new conflict nodes by checking the vv_locks
            // vv_locks stores all the vertices in the VV(v0) & VV(v1) that are not contained by n or outer_v_block
            update_conflict_nodes(vv_locks, n_id, tree);
        }
        for (iset_iter it = vv_locks.begin(); it != vv_locks.end(); it++)
        {
            omp_unset_lock(&(v_locks[*it - 1]));
        }
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
    }

    // leaf_VV vvs;
    // n.get_VV(vvs,mesh);
}

void Contraction_Simplifier::get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1, Node_V *&outer_v_block, Node_V &n, Mesh &mesh, leaf_VT &vts, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, PRT_Tree &tree)
{

    //cout<<"[NOTICE]get edge relation"<<endl;
    outer_v_block = NULL;
    /// inverted order as I only need the block indexing v1

    if (e[1] > e[0])
    {
        vt1 = Contraction_Simplifier::get_VT(e[1], n, mesh, vts, cache, tree, outer_v_block, params);
        vt0 = Contraction_Simplifier::get_VT(e[0], n, mesh, vts, cache, tree, outer_v_block, params);
    }
    else
    {

        vt0 = Contraction_Simplifier::get_VT(e[0], n, mesh, vts, cache, tree, outer_v_block, params);
        vt1 = Contraction_Simplifier::get_VT(e[1], n, mesh, vts, cache, tree, outer_v_block, params);
    }

    Contraction_Simplifier::get_ET(e, et, n, mesh, vts);
}

//A new version of get_edge_relations for parallel computing without global cache
void Contraction_Simplifier::get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1,bool& v1_is_border, bool& v2_is_border,  Node_V *&outer_v_block, Node_V &n, Mesh &mesh, leaf_VT &vts, boost::dynamic_bitset<>is_border_edge, map<int, leaf_VT> &cache, contraction_parameters &params, PRT_Tree &tree)
{

    //cout<<"[NOTICE]get edge relation"<<endl;
    outer_v_block = NULL;

    // Using local cache

    vt1 = Contraction_Simplifier::get_VT(e[1], n, mesh, vts, cache, tree, outer_v_block, params);
    vt0 = Contraction_Simplifier::get_VT(e[0], n, mesh, vts, cache, tree, outer_v_block, params);
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

    Contraction_Simplifier::get_ET(e, et, n, mesh, vts);
}

VT *Contraction_Simplifier::get_VT(int v_id, Node_V &n, Mesh &mesh, leaf_VT &vts, LRU_Cache<int, leaf_VT> &cache,
                                   PRT_Tree &tree, Node_V *&v_block, contraction_parameters &params)
{
    int local_index;
    bool debug = false;
    if (n.indexes_vertex(v_id))
    {
        if (debug)
            cout << "[get_VT] " << v_id << " -> LOCAL VERTEX " << n << endl;

        local_index = v_id - n.get_v_start();
        VT *vt = &(vts[local_index]);
        Contraction_Simplifier::clean_coboundary(*vt, mesh);
        v_block = &n;
        return vt;
    }
    else
    {
        if (debug) // if v is external
            cout << "[get_VT] " << v_id << " -> EXTERNAL VERTEX " << n << endl;

        tree.get_leaf_indexing_vertex(tree.get_root(), v_id, v_block); // the vertex should be protected by a lock in the future
        local_index = v_id - v_block->get_v_start();

        LRU_Cache<int, leaf_VT>::mapIt it_c = cache.end();
#pragma omp critical
        {
            it_c = cache.find(v_block->get_v_start()); //First check in the cache
            if (it_c == cache.end())                   //if not in the cache
            {
                if (debug)
                    cout << "    -> LEAF BLOCK OUTSIDE CACHE - REGEN " << *v_block << endl;

                leaf_VT lVT;
                v_block->get_VT(lVT, mesh);
                it_c = cache.insert(v_block->get_v_start(), lVT);
            }
            else
            {
                if (debug)
                {
                    cout << "    -> LEAF BLOCK IN CACHE - CLEAN " << *v_block << endl;
                }

                // if(debug/* || v_id ==2355*/)
                //     cout<<"num_elem_in_vt: "<<get_num_elements_in_container_of_containers((it_c->second)[local_index])<<endl;
                //                print_container_of_containers_content("VTop(2355) ",(it_c->second)[local_index]);
            }
            Contraction_Simplifier::clean_coboundary((it_c->second)[local_index], mesh);

            if (debug)
                cout << "[NOTICE]cleaned coboundary" << endl;
        }

        return &(it_c->second)[local_index];
    }
}

// New version of get_VT() without using global cache
VT *Contraction_Simplifier::get_VT(int v_id, Node_V &n, Mesh &mesh, leaf_VT &vts, map<int, leaf_VT> &cache,
                                   PRT_Tree &tree, Node_V *&v_block, contraction_parameters &params)
{
    int local_index;
    bool debug = false;
    if (n.indexes_vertex(v_id))
    {
        if (debug)
            cout << "[get_VT] " << v_id << " -> LOCAL VERTEX " << n << endl;

        local_index = v_id - n.get_v_start();
        VT *vt = &(vts[local_index]);
        Contraction_Simplifier::clean_coboundary(*vt, mesh);
        v_block = &n;
        return vt;
    }
    else
    {
        if (debug) // if v is external
            cout << "[get_VT] " << v_id << " -> EXTERNAL VERTEX " << n << endl;

        tree.get_leaf_indexing_vertex(tree.get_root(), v_id, v_block); // the vertex should be protected by a lock in the future
        local_index = v_id - v_block->get_v_start();
        //LRU_Cache<int, leaf_VT>::mapIt it_c = cache.end();
        map<int, leaf_VT>::iterator it_c = cache.end();
        VT *vt = NULL;
        it_c = cache.find(v_block->get_v_start()); //First check in the cache
        if (it_c == cache.end())                   //if not in the cache
        {
            if (debug)
                cout << "    -> LEAF BLOCK OUTSIDE CACHE - REGEN " << *v_block << endl;

            leaf_VT lVT;
            v_block->get_VT(lVT, mesh);
            cache.insert({v_block->get_v_start(), lVT});
            vt = &(cache[v_block->get_v_start()][local_index]);
        }
        else
        {
            if (debug)
            {
                cout << "    -> LEAF BLOCK IN CACHE - CLEAN " << *v_block << endl;
            }
            vt = &((it_c->second)[local_index]);
        }
        Contraction_Simplifier::clean_coboundary(*vt, mesh);

        if (debug)
            cout << "[NOTICE]cleaned coboundary" << endl;
        return vt;
    }
}

leaf_VT &Contraction_Simplifier::get_VTS(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache,
                                         PRT_Tree &tree, contraction_parameters &params)
{
    int local_index;
    bool debug = false;

    LRU_Cache<int, leaf_VT>::mapIt it_c = cache.end();
#pragma omp critical
    {
        it_c = cache.find(n.get_v_start()); //First check in the cache
        if (it_c == cache.end())            //if not in the cache
        {
            if (debug)
                cout << "    -> LEAF BLOCK OUTSIDE CACHE - REGEN " << n << endl;
            leaf_VT local_vts;
            n.get_VT(local_vts, mesh);
            it_c = cache.insert(n.get_v_start(), local_vts);
        }
        else
        {
            if (debug)
            {
                cout << "    -> LEAF BLOCK IN CACHE - CLEAN " << n << endl;
            }
        }
    }
    return (it_c->second);
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
    tree.compact_vertices_lists(tree.get_root(), mesh, surviving_vertices);
    time.stop();
    time.print_elapsed_time("[TIME] Compact tree vertices lists: ");

    //    print_container_content("surviving vertices: ",surviving_vertices);
    //    mesh.print_mesh(cout);
    //    int a; cin>>a;

    time.start();
    //    cerr<<"[MESH] compact"<<endl;
    Mesh_Updater mu;
    cout << "number of surviving vertices:" << surviving_vertices.size() << endl;
    mu.clean_vertices_array(mesh, new_v_positions, surviving_vertices);
    cout << "number of deleted triangles:" << params.get_counter() << endl;
    /// NEW: the update_and_compact procedure check internally if we have removed all the top d-simplices
    bool all_deleted = mu.update_and_clean_triangles_arrays(mesh, new_v_positions, new_t_positions, params.get_counter());
    time.stop();
    time.print_elapsed_time("[TIME] Compact and update mesh: ");

    // cerr<<"[STAT] mesh "<<endl;
    // cerr<<"------------"<<endl;
    // mesh.print_mesh_stats(cerr);
    // cerr<<"------------"<<endl;

    time.start();
    //    cerr<<"[TREE] update indices in the tree"<<endl;
    ///TODO: Check triangle intersection before updating the tree.

    tree.update_tree(tree.get_root(), new_v_positions, new_t_positions, all_deleted);
    time.stop();
    time.print_elapsed_time("[TIME] Update tree (top-simplices): ");

    //    Reindexer r;
    //    r.reorganize_index_and_mesh(tree,mesh,false);

    //cerr << "[RAM peak] for updating the mesh and the tree: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
}

void Contraction_Simplifier::find_candidate_edges_QEM(Node_V &n, Mesh &mesh, leaf_VT &vts, edge_queue &edges, contraction_parameters &params)
{

    map<ivect, coord_type> edge_map;
    ivect e;
    int t_count = 0;
    for (RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const &t_id = itPair.first;
        bool not_valid = false;
        //#pragma omp critical
        {
            if (mesh.is_triangle_removed(*t_id))
            {
                //  cout<<"triangle removed"<<endl;
                not_valid = true;
            }
        }
        if (not_valid)
            continue;
        t_count++;
        Triangle &t = mesh.get_triangle(*t_id);

        for (int i = 0; i < t.vertices_num(); i++)
        {
            t.TE(i, e);

            int new_vertex_pos = -1;
            double error = compute_error(e[0], e[1], mesh, new_vertex_pos);
            assert(new_vertex_pos != -1);

            if (n.indexes_vertex(e[1])) // e (v1,v2) is a candidate edge if at least v2 is in n
            {

                if (new_vertex_pos == 1)
                {
                    int tmp = e[1];
                    e[1] = e[0];
                    e[0] = tmp;
                }
                map<ivect, coord_type>::iterator it = edge_map.find(e);
                // cout<<e[0]<<" and "<<e[1]<<endl;
                if (it == edge_map.end())
                {

                    //  Edge e((*it)[0],(*it)[1]);
                    edge_map[e] = error;
                    //Edge edge_obj(e[0],e[1]);
                    if (error - params.get_maximum_limit() < Zero)
                    {
                        //   cout<<"["<<e[0]-1<<","<<e[1]-1<<"]  Error will be introduced: "<<error<<endl;

                        edges.push(new Geom_Edge(e, error));
                        //     cout<<"ENQUEUE"<<endl;
                    }
                }
            }
        }
    }
    cout << "**** [Number] " << edges.size() << " edges enqueued. Start simplification.****" << endl;
    cout << "number of remaining triangles:" << t_count << endl;

    // cout<<"======NEXT NODE======"<<endl;
}

void Contraction_Simplifier::compute_initial_QEM(Mesh &mesh, vector<dvect> &planes)
{

    for (int i = 1; i <= mesh.get_triangles_num(); i++)
    {
        /* faces are triangles */
        for (int j = 0; j < 3; j++)
        {
            double *a = &(planes[i - 1][0]);
            initialQuadric[mesh.get_triangle(i).TV(j)] += Matrix(a);
        }
    }
}

void Contraction_Simplifier::compute_triangle_plane(Mesh &mesh, vector<dvect> &trPl)
{
    double coords[3][3];
    for (int i = 1; i <= mesh.get_triangles_num(); i++)
    {

        for (int v = 0; v < 3; v++)
        {
            Vertex v1 = mesh.get_vertex(mesh.get_triangle(i).TV(v));
            coords[0][v] = v1.get_x();
            coords[1][v] = v1.get_y();
            coords[2][v] = v1.get_z();
        }

        double a, b, c, m;

        a = (coords[1][1] - coords[1][0]) * (coords[2][2] - coords[2][0]) - (coords[2][1] - coords[2][0]) * (coords[1][2] - coords[1][0]);

        b = (coords[2][1] - coords[2][0]) * (coords[0][2] - coords[0][0]) - (coords[0][1] - coords[0][0]) * (coords[2][2] - coords[2][0]);

        c = (coords[0][1] - coords[0][0]) * (coords[1][2] - coords[1][0]) - (coords[1][1] - coords[1][0]) * (coords[0][2] - coords[0][0]);

        double tmp = a * a + b * b + c * c;
        m = sqrt(tmp);

        a = a / m;
        b = b / m;
        c = c / m;

        trPl[i - 1][0] = round(a * 1000000) / 1000000.0;
        trPl[i - 1][1] = round(b * 1000000) / 1000000.0;
        trPl[i - 1][2] = round(c * 1000000) / 1000000.0;
        trPl[i - 1][3] = -1 * round((a * coords[0][0] + b * coords[1][0] + c * coords[2][0]) * 1000000) / 1000000.0;
    }
}

double Contraction_Simplifier::compute_error(int v1, int v2, Mesh &mesh, int &new_vertex_pos)
{
    double min_error;
    Matrix q_bar;
    Matrix q_delta;

    /* computer quadric of virtual vertex vf */
    //cout<<"v1: "<<v1<<" and v2: "<<v2<<endl;
    q_bar = initialQuadric[v1] + initialQuadric[v2];

    q_delta = Matrix(q_bar[0], q_bar[1], q_bar[2], q_bar[3],
                     q_bar[4], q_bar[5], q_bar[6], q_bar[7],
                     q_bar[8], q_bar[9], q_bar[10], q_bar[11],
                     0, 0, 0, 1);

    Vertex vertex_1 = mesh.get_vertex(v1);
    double vx1 = vertex_1.get_x();
    double vy1 = vertex_1.get_y();
    double vz1 = vertex_1.get_z();

    Vertex vertex_2 = mesh.get_vertex(v2);
    double vx2 = vertex_2.get_x();
    double vy2 = vertex_2.get_y();
    double vz2 = vertex_2.get_z();

    double error1 = vertex_error(q_bar, vx1, vy1, vz1);
    double error2 = vertex_error(q_bar, vx2, vy2, vz2);

    min_error = std::min(error1, error2);

    if (fabs(error1 - min_error) < Zero)
    {
        new_vertex_pos = 0;
        min_error = error1;
    }
    else
    {
        new_vertex_pos = 1;
        min_error = error2;
    }

    //  min_error = vertex_error(q_bar, new_vertex[0], new_vertex[1], new_vertex[2]);
    if (min_error < Zero)
        min_error = 0;
    return min_error;
}

void Contraction_Simplifier::preprocess(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli)
{
    //////TODO: add a timer here.
    conflict_leafs.assign(tree.get_leaves_number(), iset());
    map<int, ivect> nodes_of_t;
    v_in_leaf.assign(mesh.get_vertices_num(), -1);
    // Can be parallel
    for (int i = 0; i < tree.get_leaves_number(); i++)
    {
        Node_V *n = tree.get_leaf(i);

        // Visit all triangles in the current node
        // update nodes_of_t array
        for (RunIteratorPair itPair = n->make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const &t_id = itPair.first;
            if (n->partial_indexes_triangle_vertices(mesh.get_triangle(*t_id)))
            {
                nodes_of_t[*t_id].push_back(i);
            }
        }

        // Visit all vertices in the current leaf node
        // update v_in_leaf array
        for (auto it = n->v_array_begin_iterator(); it != n->v_array_end_iterator(); it++)
        {
            v_in_leaf[*it] = i;
        }
    }

    for (auto it = nodes_of_t.begin(); it != nodes_of_t.end(); it++)
    {
        for (int i = 0; i < it->second.size(); i++)
        {
            for (int j = i + 1; j < it->second.size(); j++)
            {
                conflict_leafs[it->second[i]].insert(it->second[j]);
                conflict_leafs[it->second[j]].insert(it->second[i]);
            }
        }
    }

    // cout<<"[DEBUG] print the conflict list"<<endl;
    // for(int i=0;i<conflict_leafs.size();i++){
    //     cout<<"The conflict list of node "<<i<<" :";
    //     for(iset_iter it=conflict_leafs[i].begin();it!=conflict_leafs[i].end();it++){
    //         cout<<*it<<", ";
    //     }
    //     cout<<endl;
    // }
}

void Contraction_Simplifier::update_conflict_nodes(VV &vv_locks, int n_id, PRT_Tree &tree)
{

    // Question: How to find the leaf node index of the node containing v
    for (auto it = vv_locks.begin(); it != vv_locks.end(); it++)
    {
        int indexed_node_id;
        //indexed_node_id = tree.get_leaf_indexing_vertex_index(tree.get_root(), vv_locks[i]);
        indexed_node_id = v_in_leaf[*it];
        if (indexed_node_id == -1)
        {
            cout << "SHOULD NOT HAPPEN" << endl;
            int a;
            cin >> a;
        }
        if (conflict_leafs[n_id].find(indexed_node_id) == conflict_leafs[n_id].end())
        {
            conflict_leafs[n_id].insert(indexed_node_id);
            conflict_leafs[indexed_node_id].insert(n_id);
        }
    }
}

void Contraction_Simplifier::generate_v_in_leaf(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli)
{
    itype num_leaf = tree.get_leaves_number();
    for (int i = 0; i < num_leaf; i++)
    {
        Node_V *leaf = tree.get_leaf(i);
        for (auto it = leaf->v_array_begin_iterator(); it != leaf->v_array_end_iterator(); it++)
        {
            v_in_leaf[*it] = i;
        }
    }
}

bool Contraction_Simplifier::valid_boundary_condition(int v1,int v2,VT &vt1, VT &vt2,ET& et, bool v1_is_border, bool v2_is_border, Mesh &mesh){

  if(v1_is_border||v2_is_border){
      //cout<<"border edge"<<endl;
    return false;}
    if(vt1.size()<4||vt2.size()<4){
       // cout<<"less than 4 triangles"<<endl;
        return false;
    }


return true;
}

bool Contraction_Simplifier::not_fold_over(int v1,int v2,VT &vt1, VT &vt2,ET& et, Mesh &mesh){
VT vt2_sub_et=vt2;
    ivect et_vec;
    et_vec.push_back(et.first);
    if (et.second != -1)
        et_vec.push_back(et.second);
difference_of_vectors(vt2_sub_et,et_vec);
for (int i=0; i<vt2_sub_et.size();i++){
    Triangle t= mesh.get_triangle(vt2_sub_et[i]);
    int v2_pos = t.vertex_index(v2);
    ivect e;
    t.TE(v2_pos,e);
    bool same_side = Geometry_Wrapper::same_side_of_edge(v1,v2,e[0],e[1],mesh);
    if (!same_side)
        return false;
}
return true;
}