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

#include "topological_queries.h"

void Topological_Queries::windowed_VT_Leaf(Node_T &n, Box &dom, Box &b, Mesh& mesh, wVT &vt)
{
    itype v_start;
    itype v_end;

    n.get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;

    leaf_VT local_vt;  // local smaller structure... in the end inserted into the global map..
    local_vt.assign(v_end-v_start,VT());

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            itype real_v_index = t.TV(v);
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (n.indexes_vertex(v_start,v_end,real_v_index) && b.contains_with_all_closed_faces(mesh.get_vertex(real_v_index)))
                local_vt[real_v_index-v_start].push_back(*t_id);
        }
    }

    // finally put the local VT into the global associative array
    for(unsigned i=0; i<local_vt.size(); i++)
    {
        if(local_vt[i].size() > 0)
        {
            itype real_v_index = i+v_start;
            vt.insert(make_pair(real_v_index,local_vt[i]));
        }
    }
}

void Topological_Queries::windowed_VT_Leaf(Node_V& n, Box &b, Mesh& mesh, wVT &vt)
{
    if(!n.indexes_vertices())
        return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();

    leaf_VT local_vt;  // local smaller structure... in the end inserted into the global map..
    local_vt.assign(v_end-v_start,VT());

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            itype real_v_index = t.TV(v);
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (n.indexes_vertex(real_v_index) && b.contains_with_all_closed_faces(mesh.get_vertex(real_v_index)))
                local_vt[real_v_index-v_start].push_back(*t_id);
        }
    }

    // finally put the local VT into the global associative array
    for(unsigned i=0; i<local_vt.size(); i++)
    {
        if(local_vt[i].size() > 0)
        {
            itype real_v_index = i+v_start;
            vt.insert(make_pair(real_v_index,local_vt[i]));
        }
    }
}

void Topological_Queries::update_resulting_VT(itype v, itype t, wVT &vt)
{
    wVT::iterator iter = vt.find(v);
    if(iter == vt.end())
    {
        // to insert into the map
        VT local;
        local.push_back(t);
        vt.insert(make_pair(v,local));
    }
    else
    {
        // add the triangle to the current
        iter->second.push_back(t);
    }
}

void Topological_Queries::add_edges(itype t_id, vector<edge_triangle_tuple> &tuples, Mesh &mesh, wTT::const_iterator &iter, wTT &tt)
{
    Triangle& t = mesh.get_triangle(t_id);

    edge_triangle_tuple tuple;
    if(iter==tt.end()) // no entries in tt map
    {
        for(int v=0; v<t.vertices_num(); v++)
        {
            t.edge_tuple(v,tuple,t_id);
            tuples.push_back(tuple);
        }
    }
    else
    {
        const ivect &partial_tt = iter->second;
        for(unsigned i=0; i<partial_tt.size(); i++)
        {
            if(partial_tt[i]==-1) //the adj is unset
            {
                t.edge_tuple(i,tuple,t_id);
                tuples.push_back(tuple);
            }
        }
    }
}

void Topological_Queries::pair_adjacent_triangles(vector<edge_triangle_tuple> &tuples, Mesh &, wTT &tt)
{
    unsigned j=0;
    while(j<tuples.size())
    {
        if(j+1<tuples.size())
        {
            if(tuples[j] == tuples[j+1])
            {
                update_resulting_TT(tuples[j].get_f_pos(),tuples[j].get_t(),tuples[j+1].get_t(),tt);
                update_resulting_TT(tuples[j+1].get_f_pos(),tuples[j+1].get_t(),tuples[j].get_t(),tt);
                j+=2;
            }
            else
            {
                j++;
            }
        }
        else
        {
            j++;
        }

        if(j>=tuples.size())
            break;
    }
}

void Topological_Queries::update_resulting_TT(int pos, itype t1, itype t2, wTT &tt)
{
    wTT::iterator iter = tt.find(t1);
    if(iter == tt.end())
    {
        cout<<"[update_resulting_TT] something wrong goes here..."<<endl;
        int a; cin>>a;
    }
    else
    {
        // add the triangle to the current
        iter->second[pos] = t2;
    }
}

void Topological_Queries::init_TT_entry(itype t1, utype size, wTT &tt)
{
    ivect local;
    local.assign(size,-1);
    tt.insert(make_pair(t1,local));
}

void Topological_Queries::finalize_TT_Leaf(vector<edge_triangle_tuple> &tuples, wTT &tt, Mesh &mesh)
{
    // order the faces array
    sort_edges_tuples(tuples);
    // set the adjacencies on these faces
    pair_adjacent_triangles(tuples,mesh,tt);
}

void Topological_Queries::batched_VT_leaf(Node_V &n, Box &, Mesh &mesh, bool stats, utype &max_entries)
{
    leaf_VT local_vt;
    n.get_VT(local_vt,mesh);

    if(stats)
    {
        utype entries = 0;

        for(leaf_VT::iterator it = local_vt.begin(); it != local_vt.end(); ++it)
        {
            entries += it->size();
        }

        if(max_entries < entries)
            max_entries = entries;
    }
}

void Topological_Queries::batched_VT_leaf(Node_T &n, Box &dom, Mesh &mesh, bool stats, utype &max_entries)
{
    leaf_VT local_vt;
    n.get_VT(local_vt,dom,mesh);

    if(stats)
    {
        utype entries = 0;

        for(leaf_VT::iterator it = local_vt.begin(); it != local_vt.end(); ++it)
        {
            entries += it->size();
        }

        if(max_entries < entries)
            max_entries = entries;
    }
}

void Topological_Queries::batched_VT_no_reindex_leaf(Node_V &n, Box &dom, Mesh &mesh, bool stats, utype &max_entries)
{
    if(!n.indexes_vertices())
        return; // no vertices.. skip the current leaf block

    wTT local_vt;  // local smaller structure... in the end inserted into the global map..

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (dom.contains(mesh.get_vertex(t.TV(v)),mesh.get_domain().get_max()))
                update_resulting_VT(t.TV(v),*t_id,local_vt);
        }
    }

    if(stats)
    {
        utype entries = 0;

        for(wTT::iterator it = local_vt.begin(); it != local_vt.end(); ++it)
        {
            entries += it->second.size();
        }

        if(max_entries < entries)
            max_entries = entries;
    }
}

void Topological_Queries::batched_VT_no_reindex_leaf(Node_T &n, Box &dom, Mesh &mesh, bool stats, utype &max_entries)
{
    wTT local_vt;  // local smaller structure... in the end inserted into the global map..

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (dom.contains(mesh.get_vertex(t.TV(v)),mesh.get_domain().get_max()))
                update_resulting_VT(t.TV(v),*t_id,local_vt);
        }
    }

    if(stats)
    {
        utype entries = 0;

        for(wTT::iterator it = local_vt.begin(); it != local_vt.end(); ++it)
        {
            entries += it->second.size();
        }

        if(max_entries < entries)
            max_entries = entries;
    }
}
