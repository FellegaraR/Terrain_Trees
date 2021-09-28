/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Aspect.h
 * Author: ytsong
 *
 * Created on January 23, 2019, 3:45 PM
 */

#ifndef ASPECT_H
#define ASPECT_H

#include "terrain_trees/node_v.h"
#include "terrain_trees/node_t.h"
#include "terrain_trees/spatial_subdivision.h"
class Aspect {
public:
    Aspect()
    {        min = INT_MAX;
        max = INT_MIN;
        avg = 0;
        num = 0;}
    void compute_triangles_aspects(Node_V &n, Mesh &mesh, Spatial_Subdivision &division);
    void compute_triangles_aspects(Node_T &n, Box &dom, int level, Mesh &mesh, Spatial_Subdivision &division);
    template<class N> void compute_triangles_aspects(N &n, Mesh &mesh, Spatial_Subdivision &division);
    inline void print_aspects_stats() { cerr<<"   min: "<<min<<" avg: "<<avg/(coord_type)num<<" max: "<<max<<endl<<" num_edges: "<<num<<endl;; }
    inline void print_aspects_stats(utype tnum) { cerr<<"   min: "<<min<<" avg: "<<avg/(coord_type)tnum<<" max: "<<max<<endl; }

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
    coord_type compute_triangle_aspect(Triangle &t, Mesh &mesh);
    void triangle_aspects_leaf(Node_V& n, Mesh& mesh);
    void triangle_aspects_leaf(Node_T& n, Box &dom, Mesh& mesh);
    template<class N> void triangle_aspects_leaf(N& n, Mesh& mesh);
};

template<class N> void Aspect::compute_triangles_aspects(N &n, Mesh &mesh, Spatial_Subdivision &division)
{
    if (n.is_leaf())
    {
        this->triangle_aspects_leaf(n,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->compute_triangles_aspects(*n.get_son(i),mesh,division);
            }
        }
    }
}

template<class N> void Aspect::triangle_aspects_leaf(N& n, Mesh& mesh)
{
    dvect aspects;
    aspects.assign(n.get_real_t_array_size(),0);
    int t_counter = 0;

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        aspects[t_counter] =this->compute_triangle_aspect(t,mesh);
        t_counter++;
    }
}
#endif /* ASPECT_H */

