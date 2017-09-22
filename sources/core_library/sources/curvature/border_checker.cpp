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

#include "border_checker.h"

void Border_Checker::compute_borders(Node_T &n, Box &dom, int level, Mesh &mesh, Spatial_Subdivision &division)
{
    if (n.is_leaf())
    {
        this->compute_borders(n,dom,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->compute_borders(*n.get_son(i),son_dom,son_level,mesh,division);
        }
    }
}

void Border_Checker::compute_borders(Node_V &n, Box &dom, int level, Mesh &mesh, Spatial_Subdivision &division)
{
    if (n.is_leaf())
    {
        this->compute_borders(n,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            this->compute_borders(*n.get_son(i),dom,level,mesh,division);
        }
    }
}

void Border_Checker::compute_borders(Node_V &n, Mesh &mesh)
{
    if(!n.indexes_vertices())
        return;

    vector< vector<edge_triangle_tuple> > all_edges;
    vector<edge_triangle_tuple> edges;
    all_edges.assign(n.get_v_end()-n.get_v_start(),edges);

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(itype j=0; j<t.vertices_num(); j++)
        {
            itype real_index = t.TV(j);
            if(n.indexes_vertex(real_index))
            {
                //we insert the three triangular faces incident in vertex v_pos
                get_incident_edges(t,*t_id,j,all_edges[real_index-n.get_v_start()]);
            }
        }
    }

    for(auto etuples : all_edges)
    {
        if(etuples.size() > 0)
        {
            this->set_borders(etuples,mesh);
        }
    }
}

void Border_Checker::compute_borders(Node_T &n, Box &n_dom, Mesh& mesh)
{
    itype v_start;
    itype v_end;

    n.get_v_range(v_start,v_end,n_dom,mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;

    vector< vector<edge_triangle_tuple> > all_edges;
    vector<edge_triangle_tuple> edges;
    all_edges.assign(v_end-v_start,edges);

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int j=0; j<t.vertices_num(); j++)
        {
            itype real_index = t.TV(j);
            if(n.indexes_vertex(v_start,v_end,real_index))
            {
                //we insert the three triangular faces incident in vertex v_pos
                get_incident_edges(t,*t_id,j,all_edges[real_index-v_start]);
            }
        }
    }

    for(auto etuples : all_edges)
    {
        if(etuples.size() > 0)
        {
            this->set_borders(etuples,mesh);
        }
    }
}

bool Border_Checker::set_borders(vector<edge_triangle_tuple> &edges, Mesh &mesh)
{
    bool borderChange = false;

    sort_edges_tuples(edges);

    unsigned j=0;
    while(j<edges.size())
    {
        if(j+1<edges.size())
        {
            if(edges[j] != edges[j+1])
            {
                Triangle& t = mesh.get_triangle(edges[j].get_t());

                for(int v=0;v<t.vertices_num();v++)
                {
                    itype v_ind = abs(t.TV(v));
                    if(edges[j].has_not(v_ind))
                    {
                        if(!t.is_border_edge(v)) //identified a triangular edge on the mesh border
                            borderChange = true;
                        t.setTV(v,-v_ind); //we flag with a negative index the vertex opposite to the edge on the border
                        break;
                    }
                }

                j++;

            }
            else if(edges[j] == edges[j+1])
            {
                j+=2;
            }
        }
        else
        {
            Triangle& t = mesh.get_triangle(edges[j].get_t());

            for(int v=0;v<t.vertices_num();v++)
            {
                itype v_ind = t.TV(v);
                if(edges[j].has_not(v_ind))
                {
                    if(!t.is_border_edge(v)) //identified a triangular edge on the mesh border
                        borderChange = true;
                    t.setTV(v,-v_ind);
                    break;
                }
            }

            j++;

        }

        if(j>=edges.size())
            break;
    }

    return borderChange;
}

void Border_Checker::get_incident_edges(Triangle &t, itype t_id, int v_pos, vector<edge_triangle_tuple> &edges)
{
    //we insert the three triangular faces incident in vertex v_pos
    edge_triangle_tuple new_item;

    for(int i=1; i<t.vertices_num(); i++)
    {
        new_item.sort_and_set(t.TV(v_pos),t.TV((v_pos+i)%t.vertices_num()),t_id);
        edges.push_back(new_item);
    }
}
