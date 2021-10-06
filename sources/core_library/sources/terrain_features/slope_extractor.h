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

#ifndef SLOPE_EXTRACTOR_H
#define SLOPE_EXTRACTOR_H

#include "geometry/geometry_slope.h"
#include "terrain_trees/node_v.h"
#include "terrain_trees/node_t.h"
#include "terrain_trees/spatial_subdivision.h"

class Slope_Extractor
{
public:

    Slope_Extractor()
    {
        min = INT_MAX;
        max = INT_MIN;
        avg = 0;
        num = 0;

    }

    //compute for all triangles the slope
    //storing the values in a global array
    void compute_triangles_slopes(Node_V &n, Mesh &mesh, Spatial_Subdivision &division);
    void compute_triangles_slopes(Node_T &n, Box &dom, int level, Mesh &mesh, Spatial_Subdivision &division);
    void compute_triangles_slopes_new(Node_V &n, Mesh &mesh, Spatial_Subdivision &division);
    void compute_triangles_slopes_new(Node_T &n, Box &dom, int level, Mesh &mesh, Spatial_Subdivision &division);


    // compute for all the triangles the slopes - using a local data structure within each leaf block
    template<class N> void compute_triangles_slopes(N &n, Mesh &mesh, Spatial_Subdivision &division);

        // compute for all the triangles the slopes - using a local data structure within each leaf block
    template<class N> void compute_triangles_slopes_new(N &n, Mesh &mesh, Spatial_Subdivision &division);

    //compute the slope values for all the edges
    void compute_edges_slopes(Node_V &n, Mesh &mesh, Spatial_Subdivision &division/*, bool set_avg_v_slopes*/);
    void compute_edges_slopes(Node_T &n, Box &dom, int level, Mesh &mesh, Spatial_Subdivision &division);

    inline void print_slopes_stats() { cerr<<"   min: "<<min<<" avg: "<<avg/(coord_type)num<<" max: "<<max<<endl<<" num_edges: "<<num<<endl;; }
    inline void print_slopes_stats(utype tnum) { cerr<<"   min: "<<min<<" avg: "<<avg/(coord_type)tnum<<" max: "<<max<<endl; }

    inline void reset_stats()
    {
        min = INT_MAX;
        max = INT_MIN;
        avg = 0;
        num = 0;
    }

private:
    coord_type min, avg, max;
    int num;

    void triangle_slopes_leaf(Node_V& n, Mesh& mesh);
    void triangle_slopes_leaf(Node_T& n, Box &dom, Mesh& mesh);
    template<class N> void triangle_slopes_leaf(N& n, Mesh& mesh);
    void triangle_slopes_leaf_new(Node_V& n, Mesh& mesh);
    void triangle_slopes_leaf_new(Node_T& n, Box &dom, Mesh& mesh);
    template<class N> void triangle_slopes_leaf_new(N& n, Mesh& mesh);
};

template<class N> void Slope_Extractor::compute_triangles_slopes(N &n, Mesh &mesh, Spatial_Subdivision &division)
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

template<class N> void Slope_Extractor::compute_triangles_slopes_new(N &n, Mesh &mesh, Spatial_Subdivision &division)
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


template<class N> void Slope_Extractor::triangle_slopes_leaf(N& n, Mesh& mesh)
{
    dvect slopes;
    slopes.assign(n.get_real_t_array_size(),0);
    int t_counter = 0;

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        slopes[t_counter] = Geometry_Slope::compute_triangle_slope(t,mesh);
        t_counter++;
    }
}

template<class N> void Slope_Extractor::triangle_slopes_leaf_new(N& n, Mesh& mesh)
{
    dvect slopes;
    slopes.assign(n.get_real_t_array_size(),0);
    int t_counter = 0;

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        slopes[t_counter] = Geometry_Slope::compute_triangle_slope(t,mesh);
        t_counter++;
    }
}
#endif // SLOPE_EXTRACTOR_H
