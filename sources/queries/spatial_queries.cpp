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

#include "spatial_queries.h"

bool Spatial_Queries::atomic_point_in_triangle_test(itype t_id, Point& p, QueryStatistics& qS, Mesh& mesh)
{
    qS.numGeometricTest++;
    if(Geometry_Wrapper::point_in_triangle(t_id,p,mesh))
    {
        qS.triangles.push_back(t_id);
        return true;
    }
    return false;
}

void Spatial_Queries::atomic_triangle_in_box_test(itype t_id, Box &b, QueryStatistics& qS, Mesh& mesh, bool get_stats)
{
    if(get_stats)
        qS.increase_tri_counter(t_id);

    if(!qS.is_checked(t_id))
    {
        qS.set_checked(t_id);

        if(get_stats)
            qS.numGeometricTest++;

        if (Geometry_Wrapper::triangle_in_box(t_id,b,mesh))
        {
            qS.triangles.push_back(t_id);
        }
    }
}

void Spatial_Queries::insert_nearest_vertices(priority_iNN_queue<Node_V> &queue, PQueue_Element<Node_V> &element, Point &p, Mesh &mesh)
{
    if(!element.get_node()->indexes_vertices())
        return; //no internal vertices..

    for(RunIteratorPair itPair = element.get_node()->make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& v_id = itPair.first;
        Vertex &v = mesh.get_vertex(*v_id);
        coord_type v_dist = v.distance(p);
        if(v_dist >= element.get_distance())
        {
            PQueue_Element<Node_V> addv = PQueue_Element<Node_V>(v_dist,*v_id);
            queue.push(addv);
        }
    }
}

void Spatial_Queries::insert_nearest_vertices(priority_iNN_queue<Node_T> &queue, PQueue_Element<Node_T> &element, Point &p, Mesh &mesh)
{
    itype v_start;
    itype v_end;

    element.get_node()->get_v_range(v_start,v_end,element.get_domain(),mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;

    for(itype v_id=v_start; v_id < v_end; v_id++)
    {
        Vertex &v = mesh.get_vertex(v_id);
        coord_type v_dist = v.distance(p);
        if(v_dist >= element.get_distance())
        {
            PQueue_Element<Node_T> addv = PQueue_Element<Node_T>(v_dist,v_id);
            queue.push(addv);
        }
    }
}
