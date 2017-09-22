#include "pmrt_tree.h"

PMRT_Tree::PMRT_Tree(int triangles_per_leaf, int sons_num) : Tree<Node_T>(sons_num)
{
    this->triangles_threshold = triangles_per_leaf;
    this->mesh = Mesh();
    this->root = Node_T();
}

void PMRT_Tree::build_tree()
{
    for (itype i = 1; i <= this->mesh.get_triangles_num(); i++)
    {
        this->add_triangle(this->root, this->mesh.get_domain(), 0, i);
    }
}

void PMRT_Tree::add_triangle(Node_T& n, Box& domain, int level, itype t)
{
    if (!Geometry_Wrapper::triangle_in_box_build(t,domain,this->mesh)) return;

    if (n.is_leaf())
    {
        n.add_triangle(t);
        if (is_full(n))
        {
            this->split(n, domain, level);
        }
    }
    else
    {
        for (int i = 0; i<this->subdivision.son_number(); i++)
        {
            Box son_dom = this->subdivision.compute_domain(domain, level, i);
            int son_level = level +1;
            this->add_triangle(*n.get_son(i), son_dom, son_level, t);
        }
    }
}

void PMRT_Tree::reinsert_triangle_once(Node_T& n, Box& domain, itype t)
{
    if (!Geometry_Wrapper::triangle_in_box_build(t,domain,this->mesh)) return;

    if (n.is_leaf())
    {
        n.add_triangle(t);
    }
}

void PMRT_Tree::split(Node_T& n, Box& domain, int level)
{
    n.init_sons(this->subdivision.son_number());

    for (int i = 0; i<this->subdivision.son_number(); i++)
    {
        Node_T* s = new Node_T();
        n.set_son(s,i);
        s->set_parent(&n);
    }

    // we reinsert the triangles in the son nodes only once
    // thus, without splitting any further the space
    for (int j = 0; j<this->subdivision.son_number(); j++)
    {
        Box son_dom = this->subdivision.compute_domain(domain,level, j);
        for(RunIterator runIt = n.t_array_begin_iterator(), runEnd = n.t_array_end_iterator(); runIt != runEnd; ++runIt)
        {
            this->reinsert_triangle_once(*n.get_son(j) , son_dom, *runIt);
        }
    }

    n.clear_t_array();
}
