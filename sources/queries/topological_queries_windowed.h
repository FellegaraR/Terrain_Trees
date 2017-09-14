/*
    This file is part of the Terrain Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The Terrain Trees library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Terrain Trees library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Terrain Trees library.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "topological_queries.h"
#include "utilities/timer.h"
#include "io/reader.h"
#include "utilities/sorting.h"
#include "io/writer.h"

template<class N> void Topological_Queries::windowed_VT(N &n, Box &dom, Mesh &mesh, Spatial_Subdivision &division, string query_path, bool reindexed)
{
    vector<Box> boxes;
    Reader::read_queries(boxes,query_path);

    wVT results;

    Timer time;
    coord_type tot_time = 0;


    for(unsigned j=0;j<boxes.size();j++)
    {
        time.start();
        if(reindexed)
            windowed_VT(n,dom,0,boxes[j],mesh,division,results);
        else
            windowed_VT_no_reindex(n,dom,0,boxes[j],mesh,division,results);
        time.stop();
        tot_time += time.get_elapsed_time();

        //debug print
        cout<<"for box "<<j<<" vertices found: "<<results.size()<<endl;
        results.clear();
    }    
    cerr<<"[TIME] extracting windowed VT "<<tot_time<<endl;
}

template<class D> void Topological_Queries::windowed_VT(Node_T &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, wVT &vt)
{
    if (!dom.intersects(b))
        return;

    if (n.is_leaf())
    {
        this->windowed_VT_Leaf(n,dom,b,mesh,vt);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->windowed_VT(*n.get_son(i), son_dom, son_level, b, mesh, division, vt);
        }
    }
}

template<class N, class D> void Topological_Queries::windowed_VT_no_reindex(N &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, wVT &vt)
{
    if (!dom.intersects(b))
        return;

    if (n.is_leaf())
    {
        this->windowed_VT_Leaf_no_reindex(n,dom,b,mesh,vt);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->windowed_VT_no_reindex(*n.get_son(i), son_dom, son_level, b, mesh, division, vt);
        }
    }
}

template<class D> void Topological_Queries::windowed_VT(Node_V &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, wVT &vt)
{
    if (!dom.intersects(b))
        return;

    if (n.is_leaf())
    {
        // here we do not need the son dom to understand if a vertex is inside the leaf ==>> we use the vertices range
        this->windowed_VT_Leaf(n,b,mesh,vt);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->windowed_VT(*n.get_son(i), son_dom, son_level, b, mesh, division, vt);
        }
    }
}

template<class N> void Topological_Queries::windowed_VT_Leaf_no_reindex(N& n, Box &dom, Box &b, Mesh& mesh, wVT &vt)
{
    wVT local_vt;  // local smaller structure... in the end inserted into the global map..

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (dom.contains(mesh.get_vertex(t.TV(v)),mesh.get_domain().get_max()) &&
                    b.contains_with_all_closed_faces(mesh.get_vertex(t.TV(v))))
                update_resulting_VT(t.TV(v),*t_id,local_vt);
        }
    }

    vt.insert(local_vt.begin(),local_vt.end());
}

template<class N> void Topological_Queries::windowed_TT(N &n, Box &dom, Mesh &mesh, Spatial_Subdivision &division, string query_path)
{
    vector<Box> boxes;
    Reader::read_queries(boxes,query_path);

    wTT results;
    Timer time;
    coord_type tot_time = 0.0;

    boost::dynamic_bitset<> checked_tri(mesh.get_triangles_num());

    for(unsigned j=0;j<boxes.size();j++)
    {
        time.start();
        windowed_TT(n,dom,0,boxes[j],mesh,division,results,checked_tri);
        time.stop();
        tot_time += time.get_elapsed_time();

        //debug print
        cout<<"for box "<<j<<" tetrahedra found: "<<results.size()<<endl;
        results.clear();

        checked_tri.reset();
    }
    cerr<<"[TIME] extracting windowed TT "<<tot_time<<endl;
}

template<class N, class D> void Topological_Queries::windowed_TT(N &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, wTT &tt,
                                                                 boost::dynamic_bitset<> &checked_tri)
{
    if (!dom.intersects(b))
        return;

    if (n.is_leaf())
    {
        if(b.completely_contains(dom))
            this->windowed_TT_Leaf_add(n,mesh,tt,checked_tri);
        else
            this->windowed_TT_Leaf_test(n,b,mesh,tt,checked_tri);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->windowed_TT(*n.get_son(i), son_dom, son_level, b, mesh, division, tt, checked_tri);
        }
    }
}

template<class N> void Topological_Queries::windowed_TT_Leaf_test(N& n, Box &b, Mesh& mesh, wTT &tt, boost::dynamic_bitset<> &checked_tri)
{
    vector<edge_triangle_tuple> faces;
    Box bb;
    pair<itype,itype> run;

    for(ivect::iterator it=n.get_t_array_begin(); it!=n.get_t_array_end(); ++it)
    {
        if(n.get_run_bounding_box(it,bb,mesh,run))
        {
            if(b.completely_contains(bb))
            {
                for(itype t_id=run.first; t_id<=run.second; t_id++)
                {
                    wTT::const_iterator entry = tt.find(t_id);

                    checked_tri[t_id-1]=true;

                    //if the run is completely contained.. simply add..
                    if(entry == tt.end()) // first time for the current triangle
                        init_TT_entry(t_id,mesh.get_triangle(t_id).vertices_num(),tt);
                    add_edges(t_id,faces,mesh,entry,tt);
                }
            }
            else if(b.intersects(bb))
            {
                for(int t_id=run.first; t_id<=run.second; t_id++)
                {
                    wTT::const_iterator entry = tt.find(t_id);

                    //if I have an entry into the result or I have an intersection with the box
                    if(entry != tt.end() || (!checked_tri[t_id-1] && Geometry_Wrapper::triangle_in_box(t_id,b,mesh)))
                    {

                        if(entry == tt.end()) // first time for the current triangle
                            init_TT_entry(t_id,mesh.get_triangle(t_id).vertices_num(),tt);
                        add_edges(t_id,faces,mesh,entry,tt);
                    }

                    checked_tri[t_id-1]=true;
                }
            }
        }
        else
        {
            wTT::const_iterator entry = tt.find(*it);

            //if I have an entry into the result or I have an intersection with the box
            if(entry != tt.end() || (!checked_tri[*it-1] && Geometry_Wrapper::triangle_in_box(*it,b,mesh)))
            {

                if(entry == tt.end()) // first time for the current triangle
                    init_TT_entry(*it,mesh.get_triangle(*it).vertices_num(),tt);
                add_edges(*it,faces,mesh,entry,tt);
            }

            checked_tri[*it-1]=true;
        }
    }
    finalize_TT_Leaf(faces,tt,mesh);
}

template<class N> void Topological_Queries::windowed_TT_Leaf_add(N& n, Mesh& mesh, wTT &tt, boost::dynamic_bitset<> &checked_tri)
{
    vector<edge_triangle_tuple> faces;

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;

        wTT::const_iterator entry = tt.find(*t_id);

        checked_tri[*t_id-1]=true;

        //if I have an entry into the result or I have an intersection with the box
        if(entry == tt.end()) // first time for the current triangle
            init_TT_entry(*t_id,mesh.get_triangle(*t_id).vertices_num(),tt);
        add_edges(*t_id,faces,mesh,entry,tt);
    }
    finalize_TT_Leaf(faces,tt,mesh);
}
