#include "spatial_subdivision.h"
#include <boost/dynamic_bitset.hpp>

Box Spatial_Subdivision::compute_domain(Box& parent_dom, int level, int child_ind)
{
    if(sons == 4) //quadtree
        return compute_quad_domain(parent_dom,child_ind);
    else if(sons == 2) // kD-tree
        return compute_kd_domain(parent_dom,level,child_ind);
    else
    {
        cerr<<"[ERROR] invalid spatial decomposition"<<endl;
        return Box();
    }
}

Box Spatial_Subdivision::compute_quad_domain(Box& parent_dom, int child_ind)
{
    boost::dynamic_bitset<> son_id_bits(2,child_ind); // 2D
    Box sonDom = Box();

    Point &min = sonDom.get_min();
    Point &max = sonDom.get_max();

    Point &p_min = parent_dom.get_min();
    Point &p_max = parent_dom.get_max();

    for(int i=0; i<min.get_dimension(); i++) // 2D
    {
        if(son_id_bits[i])
        {
            min.set_c(i,p_min.get_c(i)+(p_max.get_c(i)-p_min.get_c(i))/2.0);
            max.set_c(i,p_max.get_c(i));
        }
        else
        {
            min.set_c(i,p_min.get_c(i));
            max.set_c(i,p_min.get_c(i)+(p_max.get_c(i)-p_min.get_c(i))/2.0);
        }
    }
    return sonDom;
}

Box Spatial_Subdivision::compute_kd_domain(Box& parent_dom, int level, int child_ind)
{
    int coord_to_change = level % 2; //2D

    Box sonDom = parent_dom;

    if(child_ind == 1)
    {
        sonDom.get_min().set_c(coord_to_change,parent_dom.get_min().get_c(coord_to_change)+
                               (parent_dom.get_max().get_c(coord_to_change)-
                                parent_dom.get_min().get_c(coord_to_change))/2.0);
    }
    else if(child_ind == 0)
    {
        sonDom.get_max().set_c(coord_to_change,parent_dom.get_min().get_c(coord_to_change)+
                               (parent_dom.get_max().get_c(coord_to_change)-
                                parent_dom.get_min().get_c(coord_to_change))/2.0);
    }

    return sonDom;
}

