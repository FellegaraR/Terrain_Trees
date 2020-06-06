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

#include "node_t.h"

void Node_T::get_v_range(itype &v_start, itype &v_end, Box& dom, Mesh& mesh)
{
    v_start = v_end = -1;

    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& runIt = itPair.first;
        Triangle& t = mesh.get_triangle(*runIt);
        for(int v=0; v<t.vertices_num(); v++)
        {
            itype v_id = t.TV(v);
            if((v_start == -1 && v_end == -1) || (v_id < v_start || v_id >= v_end))
                //I check only the vertices outside the current range
            {
                if(dom.contains(mesh.get_vertex(v_id),mesh.get_domain().get_max()))
                {
                    if(v_start == -1 || v_start > v_id)
                        v_start = v_id;
                    if(v_end == -1 || v_end <= v_id)
                        v_end = v_id+1;// v_end is +1 because it represents the first index outside the leaf
                }
            }
        }
    }
}

void Node_T::get_VT(leaf_VT &all_vt, Box &dom, Mesh &mesh)
{
    itype v_start;
    itype v_end;

    get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;

    all_vt.assign(v_end-v_start,ivect());

    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            itype real_v_index = t.TV(v);
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (indexes_vertex(v_start,v_end,real_v_index))
                all_vt[real_v_index-v_start].push_back(*t_id);
        }
    }
}

void Node_T::get_VT(leaf_VT &all_vt, itype v_start, itype v_end, Mesh &mesh)
{
    all_vt.assign(v_end-v_start,ivect());

    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            itype real_v_index = t.TV(v);
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (indexes_vertex(v_start,v_end,real_v_index))
                all_vt[real_v_index-v_start].push_back(*t_id);
        }
    }
}

void Node_T::get_VV(leaf_VV &all_vv, Box &dom, Mesh& mesh)
{
    itype v_start;
    itype v_end;

    get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;

    all_vv.assign(v_end-v_start,iset());

    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            itype v_id = t.TV(v);
            if (indexes_vertex(v_start,v_end,v_id))
            {
                //init delle VV
                itype v_pos = v_id - v_start;
                for(int j=1; j<t.vertices_num(); j++)
                    all_vv[v_pos].insert(t.TV((v+j)%t.vertices_num()));
            }
        }
    }
}

void Node_T::get_VV(leaf_VV &all_vv, itype v_start, itype v_end, Mesh& mesh)
{
    all_vv.assign(v_end-v_start,iset());

    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            itype v_id = t.TV(v);
            if (indexes_vertex(v_start,v_end,v_id))
            {
                //init delle VV
                itype v_pos = v_id - v_start;
                for(int j=1; j<t.vertices_num(); j++)
                    all_vv[v_pos].insert(t.TV((v+j)%t.vertices_num()));
            }
        }
    }
}

void Node_T::get_VV_vector(leaf_VV_vec &all_vv, Box &dom, Mesh& mesh)
{
    itype v_start;
    itype v_end;

    get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;

    all_vv.assign(v_end-v_start,ivect());

    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            itype v_id = t.TV(v);
            if (indexes_vertex(v_start,v_end,v_id))
            {
                //init delle VV
                itype v_pos = v_id - v_start;
                for(int j=1; j<t.vertices_num(); j++){
                    if(t.is_border_edge((v+j)%t.vertices_num()))
                {
                    all_vv[v_pos].push_back(t.TV((v+t.vertices_num()-j)%t.vertices_num()));   
                }
                    all_vv[v_pos].push_back(t.TV((v+j)%t.vertices_num()));}
            }
        }
    }
}

void Node_T::get_VV_vector(leaf_VV_vec &all_vv, itype v_start, itype v_end, Mesh& mesh)
{
    all_vv.assign(v_end-v_start,ivect());

    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            itype v_id = t.TV(v);
            if (indexes_vertex(v_start,v_end,v_id))
            {
                //init delle VV
                itype v_pos = v_id - v_start;
                for(int j=1; j<t.vertices_num(); j++)
                {   if(t.is_border_edge((v+j)%t.vertices_num()))
                {
                    all_vv[v_pos].push_back(t.TV((v+t.vertices_num()-j)%t.vertices_num()));   
                }
                    all_vv[v_pos].push_back(t.TV((v+j)%t.vertices_num()));}
            }
        }
    }
}

void Node_T::get_VV_VT(leaf_VV &all_vv, leaf_VT &all_vt, itype v_start, itype v_end, Mesh &mesh)
{
    all_vv.assign(v_end-v_start,iset());
    all_vt.assign(v_end-v_start,ivect());

    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            itype v_id = t.TV(v);
            if (indexes_vertex(v_start,v_end,v_id))
            {
                //init delle VV
                itype v_pos = v_id - v_start;
                all_vt[v_pos].push_back(*t_id);
                for(int j=1; j<t.vertices_num(); j++)
                    all_vv[v_pos].insert(t.TV((v+j)%t.vertices_num()));
            }
        }
    }
}
