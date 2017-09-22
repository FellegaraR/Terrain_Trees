#include "abstract_curvature.h"

void Abstract_Curvature::compute(Node_V &n, Mesh &mesh, Spatial_Subdivision &division)
{
    if (n.is_leaf())
    {
        this->curvature_leaf(n,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->compute(*n.get_son(i),mesh,division);
            }
        }
    }
}

void Abstract_Curvature::compute(Node_T &n, Box &n_dom, int level, Mesh &mesh, Spatial_Subdivision &division)
{
    if (n.is_leaf())
    {
        this->curvature_leaf(n,n_dom,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(n_dom,level,i);
            int son_level = level +1;
            if(n.get_son(i)!=NULL)
            {
                this->compute(*n.get_son(i),son_dom,son_level,mesh,division);
            }
        }
    }
}
