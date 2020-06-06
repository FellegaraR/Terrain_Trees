/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Aspect.cpp
 * Author: ytsong
 * 
 * Created on January 23, 2019, 3:45 PM
 */

#include "Aspect.h"
coord_type Aspect::compute_triangle_aspect(Triangle &t, Mesh &mesh)
{
    Vertex &v1 = mesh.get_vertex(t.TV(0));
    Vertex &v2 = mesh.get_vertex(t.TV(1));
    Vertex &v3 = mesh.get_vertex(t.TV(2));

    // get the two vectors in the triangle
    dvect u = { v2.get_x()-v1.get_x() , v2.get_y()-v1.get_y() , v2.get_z()-v1.get_z() };
    dvect v = { v3.get_x()-v1.get_x() , v3.get_y()-v1.get_y() , v3.get_z()-v1.get_z() };

    // get the cross product to get the normal
    dvect n = { u[1]*v[2] - u[2]*v[1] , u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0] };
    coord_type angle;
    if (n[0]==0)
       angle=(n[1]>0)?360:180;   
    else 
    {
        angle = atan(n[1]/n[0])*180/3.14159265358979323846;
        if(n[0]>0)
            angle= 90-angle;
        else
            angle=270-angle;
    }

    

    return angle;
}

void Aspect::compute_triangles_aspects(Node_T &n, Box &dom, int level, Mesh &mesh, Spatial_Subdivision &division)
{
    if (n.is_leaf())
    {
        this->triangle_aspects_leaf(n,dom,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->compute_triangles_aspects(*n.get_son(i), son_dom, son_level, mesh, division);
        }
    }
}

void Aspect::triangle_aspects_leaf(Node_T &n, Box &dom, Mesh &mesh)
{
    itype v_start;
    itype v_end;

    n.get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;

    map<itype,coord_type> aspects;

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);

        if(n.indexes_vertex(v_start,v_end,t.minindex()))
        {
            coord_type s = this->compute_triangle_aspect(t,mesh);
            aspects[*t_id] = s;

            if(s < this->min)
                this->min = s;
            if(s > this->max)
                this->max = s;
            this->avg += s;
        }
    }
}

void Aspect::compute_triangles_aspects(Node_V &n, Mesh &mesh, Spatial_Subdivision &division)
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
void Aspect::triangle_aspects_leaf(Node_V &n, Mesh &mesh)
{
    /// if there are no vertices in the leaf we have nothing to do..
    if(!n.indexes_vertices())
        return;

    map<itype,coord_type> aspects;

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);

        if(n.indexes_vertex(t.minindex()))
        {
            coord_type s = this->compute_triangle_aspect(t,mesh);
            aspects[*t_id] = s;

            if(s < this->min)
                this->min = s;
            if(s > this->max)
                this->max = s;
            this->avg += s;
        }
    }
}