#include "node_v.h"

void Node_V::get_VT(leaf_VT &all_vt, Mesh &mesh)
{
    // here we have a reindexed index thus, if there are no vertices indexed the array size is zero
    if(!indexes_vertices())
        return;

    itype v_start = get_v_start();
    itype v_end = get_v_end();

    all_vt.assign(v_end-v_start,ivect());

    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            itype real_v_index = t.TV(v);
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (indexes_vertex(real_v_index))
                all_vt[real_v_index-v_start].push_back(*t_id);
        }
    }
}

void Node_V::get_VV(leaf_VV &all_vv,  Mesh& mesh)
{
    if(!this->indexes_vertices())
        return;

    all_vv.assign((this->get_v_end()-this->get_v_start()),iset());

    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            itype v_id = t.TV(v);
            if(this->indexes_vertex(v_id))
            {
                //init VV
                itype v_pos = v_id - this->get_v_start();
                for(int j=1; j<t.vertices_num(); j++)
                    all_vv[v_pos].insert(t.TV((v+j)%t.vertices_num()));
            }
        }
    }
}

void Node_V::get_VV_VT(leaf_VV &all_vv, leaf_VT &all_vt, Mesh &mesh)
{
    if(!this->indexes_vertices())
        return;

    itype v_start = get_v_start();
    itype v_end = get_v_end();
    all_vt.assign(v_end-v_start,ivect());
    all_vv.assign(v_end-v_start,iset());

    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            itype v_id = t.TV(v);
            if(this->indexes_vertex(v_id))
            {
                //init VV and VT
                itype v_pos = v_id - v_start;
                all_vt[v_pos].push_back(*t_id);
                for(int j=1; j<t.vertices_num(); j++)
                {
                    all_vv[v_pos].insert(t.TV((v+j)%t.vertices_num()));                    
                }

            }
        }
    }
}

void Node_V::get_ET(leaf_ET &ets,  Mesh &mesh)
{
    if(!this->indexes_vertices())
        return;

    ivect e;

    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            //ET
            t.TE(v,e);

            // an edge is considered processed only if the current leaf block indexes the extreme with the higher position index
            // in this way each edge is processed once during each traversal of the tree
            if(indexes_vertex(e[1]))
            {
                leaf_ET::iterator it = ets.find(e);
                if(it != ets.end())
                {
                    pair<itype,itype> &inside = it->second;
                    inside.second = *t_id;

                    /// force the ordering <max,min>
                    if(inside.first < inside.second)
                    {
                        itype tmp = inside.first;
                        inside.first = inside.second;
                        inside.second = tmp;
                    }
                }
                else
                {
                    /// we add to the local structure only the edges with an internal vertex
                    pair<itype,itype> new_entry = make_pair(*t_id,-1);
                    ets.insert(make_pair(e,new_entry));
                }
            }
        }
    }
}
