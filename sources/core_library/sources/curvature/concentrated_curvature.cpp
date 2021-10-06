/*
    This file is part of the Terrain Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)
               Paola Magillo (paola.magillo@unige.it)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The triangle Trees library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The triangle Trees library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the triangle Trees library.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "concentrated_curvature.h"

void Concentrated_Curvature::curvature_leaf(Node_V &n, Mesh &mesh)
{
    /// if there are no vertices in the leaf we have nothing to do..
    if(!n.indexes_vertices())
        return;

    //leaf_VT vts(n.get_v_end()-n.get_v_start(),ivect());
    boost::dynamic_bitset<> is_v_border(n.get_v_end()-n.get_v_start());
    dvect local_curvatures(n.get_v_end()-n.get_v_start(),0.0);
    vts.assign(n.get_v_end()-n.get_v_start(),ivect());
    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int j=0; j<t.vertices_num(); j++)
        {
            itype real_index = t.TV(j);
            if(n.indexes_vertex(real_index))
            {
                itype local_index = real_index-n.get_v_start();
                //popolo la VT per un dato vertice
                vts[local_index].push_back(*t_id);

                this->update_local_curvature(real_index,local_index,t,j,is_v_border,local_curvatures,mesh);
            }
        }
    }

    this->finalize_local_curvatures(n.get_v_start(),vts,is_v_border,local_curvatures,mesh);
    vts.clear();
}

void Concentrated_Curvature::curvature_leaf(Node_T &n, Box &n_dom, Mesh &mesh)
{
    itype v_start;
    itype v_end;

    n.get_v_range(v_start,v_end,n_dom,mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;

     vts.assign(v_end-v_start,ivect());
    boost::dynamic_bitset<> is_v_border(v_end-v_start);
    dvect local_curvatures(v_end-v_start,0.0);

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int j=0; j<t.vertices_num(); j++)
        {
            itype real_index = t.TV(j);
            if(n.indexes_vertex(v_start,v_end,real_index))
            {
                int local_index = real_index-v_start;
                //VT
                vts[local_index].push_back(*t_id);

                this->update_local_curvature(real_index,local_index,t,j,is_v_border,local_curvatures,mesh);
            }
        }
    }

    this->finalize_local_curvatures(v_start,vts,is_v_border,local_curvatures,mesh);

    vts.clear();
}

void Concentrated_Curvature::update_local_curvature(itype v_id, itype local_v_id, Triangle &t, itype v_pos,
                                                    boost::dynamic_bitset<> &is_v_border, dvect &local_curvatures, Mesh &mesh)
{
    itype v1 = t.TV((v_pos+1)%t.vertices_num());
    itype v2 = t.TV((v_pos+2)%t.vertices_num());
    local_curvatures[local_v_id] += acos(Geometry_Curvature::cos_angle(mesh.get_vertex(v1),mesh.get_vertex(v_id),mesh.get_vertex(v2)));

    //a vertex is on the border if at least one of the other two vertices in the triangle has a negative index
    if(!is_v_border[local_v_id])
    {
        for(int v=1; v<t.vertices_num(); v++)
        {
            if(t.is_border_edge((v_pos+v)%t.vertices_num()))
            {
                is_v_border[local_v_id] = true;
                break;
            }
        }
    }
}


void Concentrated_Curvature::finalize_local_curvatures(itype v_start, leaf_VT &vts, boost::dynamic_bitset<> &is_v_border,
                                                       dvect &local_curvatures, Mesh &mesh)
{
    for(utype i=0; i<is_v_border.size(); i++)
    {
        if (is_v_border[i])
            local_curvatures[i] = PI-local_curvatures[i]; // boundary vertex
        else
            local_curvatures[i] = 2.0*PI-local_curvatures[i];   // interior vertex

        if (type==GAUSS_ANGLE_DEFICIT)
            local_curvatures[i] /= Geometry_Curvature::fan_area(vts[i],mesh);

        mesh.get_vertex(i+v_start).add_field(local_curvatures[i]);
    }
}
