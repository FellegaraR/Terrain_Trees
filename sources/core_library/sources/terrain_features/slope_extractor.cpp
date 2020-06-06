/*
    This file is part of the Terrain Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

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
    along with the Terrain Trees library.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "slope_extractor.h"

void Slope_Extractor::compute_triangles_slopes(Node_T &n, Box &dom, int level, Mesh &mesh, Spatial_Subdivision &division)
{
    if (n.is_leaf())
    {
        this->triangle_slopes_leaf(n,dom,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->compute_triangles_slopes(*n.get_son(i), son_dom, son_level, mesh, division);
        }
    }
}

void Slope_Extractor::compute_triangles_slopes_new(Node_T &n, Box &dom, int level, Mesh &mesh, Spatial_Subdivision &division)
{
    if (n.is_leaf())
    {
        this->triangle_slopes_leaf_new(n,dom,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->compute_triangles_slopes_new(*n.get_son(i), son_dom, son_level, mesh, division);
        }
    }
}


void Slope_Extractor::triangle_slopes_leaf(Node_T &n, Box &dom, Mesh &mesh)
{
    itype v_start;
    itype v_end;

    n.get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;

    map<itype,coord_type> slopes;

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);

        if(n.indexes_vertex(v_start,v_end,t.minindex()))
        {
            coord_type s = Geometry_Slope::compute_triangle_slope(t,mesh);
            slopes[*t_id] = s;

            if(s < this->min)
                this->min = s;
            if(s > this->max)
                this->max = s;
            this->avg += s;
        }
    }
}

void Slope_Extractor::triangle_slopes_leaf_new(Node_T &n, Box &dom, Mesh &mesh)
{
    itype v_start;
    itype v_end;

    n.get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;

    map<itype,coord_type> slopes;

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);

        if(n.indexes_vertex(v_start,v_end,t.minindex()))
        {
            coord_type s = Geometry_Slope::compute_triangle_slope_song(t,mesh);
            slopes[*t_id] = s;

            if(s < this->min)
                this->min = s;
            if(s > this->max)
                this->max = s;
            this->avg += s;
        }
    }
}

void Slope_Extractor::compute_triangles_slopes(Node_V &n, Mesh &mesh, Spatial_Subdivision &division)
{
    if (n.is_leaf())
    {
        this->triangle_slopes_leaf(n,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->compute_triangles_slopes(*n.get_son(i),mesh,division);
            }
        }
    }
}

void Slope_Extractor::compute_triangles_slopes_new(Node_V &n, Mesh &mesh, Spatial_Subdivision &division)
{
    if (n.is_leaf())
    {
        this->triangle_slopes_leaf_new(n,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->compute_triangles_slopes_new(*n.get_son(i),mesh,division);
            }
        }
    }
}

void Slope_Extractor::triangle_slopes_leaf(Node_V &n, Mesh &mesh)
{
    /// if there are no vertices in the leaf we have nothing to do..
    if(!n.indexes_vertices())
        return;

    map<itype,coord_type> slopes;

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);

        if(n.indexes_vertex(t.minindex()))
        {
            coord_type s = Geometry_Slope::compute_triangle_slope(t,mesh);
            slopes[*t_id] = s;

            if(s < this->min)
                this->min = s;
            if(s > this->max)
                this->max = s;
            this->avg += s;
        }
    }
}

void Slope_Extractor::triangle_slopes_leaf_new(Node_V &n, Mesh &mesh)
{
    /// if there are no vertices in the leaf we have nothing to do..
    if(!n.indexes_vertices())
        return;

    map<itype,coord_type> slopes;

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);

        if(n.indexes_vertex(t.minindex()))
        {
            coord_type s = Geometry_Slope::compute_triangle_slope_song(t,mesh);
            slopes[*t_id] = s;

            if(s < this->min)
                this->min = s;
            if(s > this->max)
                this->max = s;
            this->avg += s;
        }
    }
}

void Slope_Extractor::compute_edges_slopes(Node_V &n, Mesh &mesh, Spatial_Subdivision &division)
{
    if (n.is_leaf())
    {
        /// if there are no vertices in the leaf we have nothing to do..
        if(!n.indexes_vertices())
            return;

        ivect e;
        map<ivect,coord_type> slopes;

        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;
            Triangle& t = mesh.get_triangle(*t_id);
            for(int i=0; i<t.vertices_num(); i++)
            {
                t.TE(i,e);
                if(n.indexes_vertex(e[0]))
                {
                    map<ivect,coord_type>::iterator it = slopes.find(e);
                    if(it == slopes.end())
                    {
                        coord_type s = Geometry_Slope::compute_edge_slope(e,mesh);
                        slopes[e] = s;

                        if(s < this->min)
                            this->min = s;
                        if(s > this->max)
                            this->max = s;
                        this->avg += s;
                        this->num++;
                    }
                }

            }
        }
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->compute_edges_slopes(*n.get_son(i),mesh,division);
            }
        }
    }
}

void Slope_Extractor::compute_edges_slopes(Node_T &n, Box &dom, int level, Mesh &mesh, Spatial_Subdivision &division)
{
    if (n.is_leaf())
    {
        itype v_start;
        itype v_end;

        n.get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range..

        if(v_start == v_end) //no internal vertices..
            return;

        ivect e;
        map<ivect,coord_type> slopes;

        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;
            Triangle& t = mesh.get_triangle(*t_id);
            for(int i=0; i<t.vertices_num(); i++)
            {
                t.TE(i,e);
                if(n.indexes_vertex(v_start,v_end,e[0]))
                {
                    map<ivect,coord_type>::iterator it = slopes.find(e);
                    if(it == slopes.end())
                    {
                        coord_type s = Geometry_Slope::compute_edge_slope(e,mesh);
                        slopes[e] = s;

                        if(s < this->min)
                            this->min = s;
                        if(s > this->max)
                            this->max = s;
                        this->avg += s;
                        this->num++;
                    }
                }

            }
        }
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->compute_edges_slopes(*n.get_son(i), son_dom, son_level, mesh, division);
        }
    }
}
