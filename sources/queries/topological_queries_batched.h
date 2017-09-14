#ifndef TOPOLOGICAL_QUERIES_BATCHED
#define TOPOLOGICAL_QUERIES_BATCHED

#include "topological_queries.h"

template<class N> void Topological_Queries::batched_VT(N &n, Box &dom, Mesh &mesh, Spatial_Subdivision &division, bool reindex)
{
    Timer time;
    utype max_entities = 0;

    time.start();
    if(reindex)
        this->batched_VT_visit(n,dom,0,mesh,division,false,max_entities);
    else
        this->batched_VT_no_reindex(n,dom,0,mesh,division,false,max_entities);
    time.stop();
    time.print_elapsed_time("[TIME] extracting bactched VT: ");

    if(reindex)
        this->batched_VT_visit(n,dom,0,mesh,division,true,max_entities);
    else
        this->batched_VT_no_reindex(n,dom,0,mesh,division,true,max_entities);
    cerr<<"[STATS] maximum number of entities: "<<max_entities<<endl;
}

template<class N, class D> void Topological_Queries::batched_VT_visit(N &n, Box &dom, int level, Mesh &mesh, D &division, bool stats, utype &max_entries)
{
    if (n.is_leaf())
    {
        this->batched_VT_leaf(n,dom,mesh,stats,max_entries);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            N& son = *n.get_son(i);
            this->batched_VT_visit(son, son_dom, son_level, mesh, division, stats, max_entries);
        }
    }
}

template<class N, class D> void Topological_Queries::batched_VT_no_reindex(N &n, Box &dom, int level, Mesh &mesh, D &division, bool stats, utype &max_entries)
{
    if (n.is_leaf())
    {
        this->batched_VT_no_reindex_leaf(n,dom,mesh,stats,max_entries);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->batched_VT_no_reindex(*n.get_son(i), son_dom, son_level, mesh, division, stats, max_entries);
        }
    }
}

template<class N> void Topological_Queries::batched_TT(N &n, Mesh &mesh, Spatial_Subdivision &division)
{
    utype max_entities = 0;

    vector<ivect> tt;
    ivect tmp;
    tmp.assign(3,-1);
    tt.assign(mesh.get_triangles_num(),tmp);

    Timer time;
    time.start();
    this->batched_TT_visit(n,mesh,division,tt,false,max_entities);
    time.stop();
    time.print_elapsed_time("[TIME] extracting bactched TT: ");

    tt.clear();
    tmp.assign(3,-1);
    tt.assign(mesh.get_triangles_num(),tmp);

    this->batched_TT_visit(n,mesh,division,tt,true,max_entities);
    cerr<<"[STATS] maximum number of faces: "<<max_entities<<endl;
}

template<class N, class D> void Topological_Queries::batched_TT_visit(N &n, Mesh &mesh, D &division, vector<ivect > &tt, bool stats, utype &max_entries)
{
    if (n.is_leaf())
    {
        this->batched_TT_leaf(n,mesh,tt,stats,max_entries);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            this->batched_TT_visit(*n.get_son(i), mesh, division, tt, stats, max_entries);
        }
    }
}

template<class N> void Topological_Queries::batched_TT_leaf(N &n, Mesh &mesh, vector<ivect > &tt, bool stats, utype &max_entries)
{
    vector<edge_triangle_tuple> faces;
    edge_triangle_tuple face;

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);

        for(int v=0; v<t.vertices_num(); v++)
        {
            if(tt[*t_id-1][v]==-1) // the entry is not initialized
            {
                t.edge_tuple(v,face,*t_id);
                faces.push_back(face);
            }
        }
    }

    // order the faces array
    sort_edges_tuples(faces);
    // set the adjacencies on these faces
    unsigned j=0;
    while(j<faces.size())
    {
        if(j+1<faces.size())
        {
            if(faces[j] == faces[j+1])
            {
                tt[faces[j].get_t() -1][faces[j].get_f_pos()] = faces[j+1].get_t();
                tt[faces[j+1].get_t() -1][faces[j+1].get_f_pos()] = faces[j].get_t();
                j+=2;
            }
            else
            {
                j++;
            }
        }
        else
        {
            j++;
        }

        if(j>=faces.size())
            break;
    }

    if(stats)
    {
        if(max_entries < faces.size())
            max_entries = faces.size();
    }
}

#endif // TOPOLOGICAL_QUERIES_BATCHED

