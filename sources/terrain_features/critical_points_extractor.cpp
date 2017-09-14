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

#include "critical_points_extractor.h"
#include "utilities/timer.h"
//#include <fstream>

void Critical_Points_Extractor::compute_critical_points(Node_V &n, Mesh &mesh, Spatial_Subdivision &division)
{
    this->critical_points.assign(mesh.get_vertices_num(),Point_Type::REGULAR); //the default is a regular vertex
    this->extract_local_critical_points(n,mesh,division);
}

void Critical_Points_Extractor::compute_critical_points(Node_T &n, Box &dom, Mesh &mesh, Spatial_Subdivision &division)
{
    this->critical_points.assign(mesh.get_vertices_num(),Point_Type::REGULAR); //the default is a regular vertex
    this->extract_local_critical_points(n,dom,0,mesh,division);
}

void Critical_Points_Extractor::extract_local_critical_points(Node_V &n, Mesh &mesh, Spatial_Subdivision &division, flat_areas &fa)
{
    if (n.is_leaf())
    {
        /// if there are no vertices in the leaf we have nothing to do..
        if(!n.indexes_vertices())
            return;

        leaf_VV vvs;
        leaf_VT vts;
        n.get_VV_VT(vvs,vts,mesh);

        this->analyze_leaf(vvs,vts,n.get_v_start(),mesh,fa);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->extract_local_critical_points(*n.get_son(i),mesh,division,fa);
            }
        }
    }
}

void Critical_Points_Extractor::extract_local_critical_points(Node_V &n, Mesh &mesh, Spatial_Subdivision &division)
{
    if (n.is_leaf())
    {
        /// if there are no vertices in the leaf we have nothing to do..
        if(!n.indexes_vertices())
            return;

        leaf_VV vvs;
        leaf_VT vts;
        n.get_VV_VT(vvs,vts,mesh);

        this->analyze_leaf(vvs,vts,n.get_v_start(),mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->extract_local_critical_points(*n.get_son(i),mesh,division);
            }
        }
    }
}

void Critical_Points_Extractor::extract_local_critical_points(Node_T &n, Box &dom, int level, Mesh &mesh, Spatial_Subdivision &division, flat_areas &fa)
{
    if(n.is_leaf())
    {
        itype v_start;
        itype v_end;

        n.get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range..

        if(v_start == v_end) //no internal vertices..
            return;

        leaf_VV vvs;
        leaf_VT vts;
        n.get_VV_VT(vvs,vts,v_start,v_end,mesh);

        this->analyze_leaf(vvs,vts,v_start,mesh,fa);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            if(n.get_son(i)!=NULL)
            {
                this->extract_local_critical_points(*n.get_son(i),son_dom,son_level,mesh,division,fa);
            }
        }
    }
}

void Critical_Points_Extractor::extract_local_critical_points(Node_T &n, Box &dom, int level, Mesh &mesh, Spatial_Subdivision &division)
{
    if(n.is_leaf())
    {
        itype v_start;
        itype v_end;

        n.get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range..

        if(v_start == v_end) //no internal vertices..
            return;

        leaf_VV vvs;
        leaf_VT vts;
        n.get_VV_VT(vvs,vts,v_start,v_end,mesh);

        this->analyze_leaf(vvs,vts,v_start,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            if(n.get_son(i)!=NULL)
            {
                this->extract_local_critical_points(*n.get_son(i),son_dom,son_level,mesh,division);
            }
        }
    }
}

void Critical_Points_Extractor::analyze_leaf(leaf_VV &vvs, leaf_VT &vts, itype v_start, Mesh &mesh, flat_areas &fa)
{
    for(unsigned i=0; i<vvs.size(); i++)
    {
        itype real_v_id = v_start + i;
        Vertex &v = mesh.get_vertex(real_v_id);
        VV &vv = vvs[i];
        bool into_flat_area = false;
        map<itype,itype> upper, lower;
        for(auto vid : vv)
        {
            Vertex &v2 = mesh.get_vertex(vid);
            if(v.get_z() == v2.get_z())
            {
                flat_areas::iterator ret = fa.find(v.get_z());
                if(ret == fa.end())
                {
                    iset triangles;
                    triangles.insert(vts[i].begin(),vts[i].end());
                    fa.insert(make_pair(v.get_z(),triangles));
                }
                else
                {
                    ret->second.insert(vts[i].begin(),vts[i].end());
                }
                into_flat_area = true;
                break;
            }
            else if(v.get_z() > v2.get_z())
            {
                lower.insert(make_pair(vid,vid)); //initially the component is formed only by the vertex itself
            }
            else
            {
                upper.insert(make_pair(vid,vid)); //initially the component is formed only by the vertex itself
            }
        }

        if(!into_flat_area)
        {
            if(upper.empty()) //the vertex is a maximum
                this->critical_points[real_v_id-1] = Point_Type::MAXIMUM;
            else if(lower.empty()) //the vertex is a minimum
                this->critical_points[real_v_id-1] = Point_Type::MINIMUM;
            else //we have to check the connected components of the edges in the link
            {
                map<itype,iset> adj_upper;
                map<itype,iset> adj_lower;

                // get the adjacent vertices that have greater or lower field values w.r.t. v
                init_adjacent_vertices(real_v_id,v,vts[i],mesh,adj_upper,adj_lower);

                // merge upper components
                utype uc_num = get_components_num(adj_upper,upper);
                // merge lower components
                utype lc_num = get_components_num(adj_lower,lower);

                if(uc_num == 2 && lc_num == 2) // simple saddle
                    this->critical_points[real_v_id-1] = Point_Type::SADDLE;
                else if(uc_num > 2 && lc_num > 2) // multiple saddle
                    this->critical_points[real_v_id-1] = Point_Type::MULTIPLE_SADDLE;
                // else the point is already flagged as regular point
            }
        }
    }
}

void Critical_Points_Extractor::analyze_leaf(leaf_VV &vvs, leaf_VT &vts, itype v_start, Mesh &mesh)
{
    for(unsigned i=0; i<vvs.size(); i++)
    {
        itype real_v_id = v_start + i;
        Vertex &v = mesh.get_vertex(real_v_id);
        VV &vv = vvs[i];
        bool into_flat_area = false;
        map<itype,itype> upper, lower;
        for(auto vid : vv)
        {
            Vertex &v2 = mesh.get_vertex(vid);
            if(v.get_z() == v2.get_z())
            {
                into_flat_area = true;
                break;
            }
            else if(v.get_z() > v2.get_z())
            {
                lower.insert(make_pair(vid,vid)); //initially the component is formed only by the vertex itself
            }
            else
            {
                upper.insert(make_pair(vid,vid)); //initially the component is formed only by the vertex itself
            }
        }

        if(!into_flat_area)
        {
            if(upper.empty()) //the vertex is a maximum
                this->critical_points[real_v_id-1] = Point_Type::MAXIMUM;
            else if(lower.empty()) //the vertex is a minimum
                this->critical_points[real_v_id-1] = Point_Type::MINIMUM;
            else //we have to check the connected components of the edges in the link
            {
                map<itype,iset> adj_upper;
                map<itype,iset> adj_lower;

                // get the adjacent vertices that have greater or lower field values w.r.t. v
                init_adjacent_vertices(real_v_id,v,vts[i],mesh,adj_upper,adj_lower);

                // merge upper components
                utype uc_num = get_components_num(adj_upper,upper);
                // merge lower components
                utype lc_num = get_components_num(adj_lower,lower);

                if(uc_num == 2 && lc_num == 2) // simple saddle
                    this->critical_points[real_v_id-1] = Point_Type::SADDLE;
                else if(uc_num > 2 && lc_num > 2) // multiple saddle
                    this->critical_points[real_v_id-1] = Point_Type::MULTIPLE_SADDLE;
                // else the point is already flagged as regular point
            }
        }
    }
}


void Critical_Points_Extractor::extract_critical_points_from_flat_areas(flat_areas &fa, Mesh &mesh)
{
    map<itype,itype> upper, lower;
    map<itype,iset> adj_upper;
    map<itype,iset> adj_lower;

    for(flat_areas::iterator it=fa.begin(); it!=fa.end(); ++it)
    {
        itype key_vertex = -1;
        coord_type field = it->first;
        iset &vt = it->second;
        for(auto tid : vt)
        {
            Triangle &t = mesh.get_triangle(tid);
            ivect link_v;
            for(int i=0; i<t.vertices_num(); i++)
            {
                Vertex &v = mesh.get_vertex(t.TV(i));
                if(v.get_z() != field)
                {
                    if(v.get_z() < field)
                        lower[t.TV(i)] = t.TV(i);
                    else
                        upper[t.TV(i)] = t.TV(i);
                    link_v.push_back(t.TV(i));
                }
                else if(key_vertex == -1)//same field and still unitialized critical vertex
                    key_vertex = t.TV(i);
            }
            if(link_v.size() == 2)
            {
                if(field < mesh.get_vertex(link_v[0]).get_z() &&
                        field < mesh.get_vertex(link_v[1]).get_z())
                {
                    adj_upper[link_v[0]].insert(link_v[1]);
                    adj_upper[link_v[1]].insert(link_v[0]);
                }
                else if(field > mesh.get_vertex(link_v[0]).get_z() &&
                        field > mesh.get_vertex(link_v[1]).get_z())
                {
                    adj_lower[link_v[0]].insert(link_v[1]);
                    adj_lower[link_v[1]].insert(link_v[0]);
                }
            }
        }

        if(upper.empty()) //the vertex is a maximum
            this->critical_points[key_vertex-1] = Point_Type::MAXIMUM;
        else if(lower.empty()) //the vertex is a minimum
            this->critical_points[key_vertex-1] = Point_Type::MINIMUM;
        else //we have to check the connected components of the edges in the link
        {
            // merge upper components
            utype uc_num = get_components_num(adj_upper,upper);
            // merge lower components
            utype lc_num = get_components_num(adj_lower,lower);

            if(uc_num == 2 && lc_num == 2) // simple saddle
                this->critical_points[key_vertex-1] = Point_Type::SADDLE;
            else if(uc_num > 2 && lc_num > 2) // multiple saddle
                this->critical_points[key_vertex-1] = Point_Type::MULTIPLE_SADDLE;
            // else the point is already flagged as regular point
        }

        // reset the local d.s.
        upper.clear();
        lower.clear();
        adj_lower.clear();
        adj_upper.clear();
    }
}

void Critical_Points_Extractor::init_adjacent_vertices(itype v_id, Vertex &v, VT &vt, Mesh &mesh,
                                                       map<itype,iset> &adj_upper, map<itype,iset> &adj_lower)
{
    ivect e;
    for(auto tid : vt)
    {
        Triangle &t = mesh.get_triangle(tid);
        int v_pos = t.vertex_index(v_id);
        t.TE(v_pos,e);

        // if the other two vertices have the fields higher/lower than v
        // -> then we flag them as adjacent
        if(v.get_z() < mesh.get_vertex(e[0]).get_z() &&
                v.get_z() < mesh.get_vertex(e[1]).get_z())
        {
            adj_upper[e[0]].insert(e[1]);
            adj_upper[e[1]].insert(e[0]);
        }
        else if(v.get_z() > mesh.get_vertex(e[0]).get_z() &&
                v.get_z() > mesh.get_vertex(e[1]).get_z())
        {
            adj_lower[e[0]].insert(e[1]);
            adj_lower[e[1]].insert(e[0]);
        }
    }
}

utype Critical_Points_Extractor::get_components_num(map<itype,iset> &adj_map, map<itype,itype> &v_flag)
{
    for(auto am : adj_map)
    {
        if(am.first == v_flag[am.first]) // not set yet
        {
            itype flag = v_flag[am.first];
            iqueue q;
            q.push(am.first);

            while(!q.empty())
            {
                itype current = q.front();
                q.pop();
                if(v_flag[current] == current) // not set
                {
                    v_flag[current] = flag;
                    iset &adj = adj_map[current];
                    for(auto a : adj)
                    {
                        q.push(a);
                    }
                }
            }
        }
    }

    iset component;
    for(map<itype,itype>::iterator it=v_flag.begin(); it!=v_flag.end(); ++it)
        component.insert(it->second);

    return component.size();
}

void Critical_Points_Extractor::print_stats()
{
    int num_reg = 0, num_min = 0, num_max = 0, num_saddle = 0, num_multisaddle = 0;
    for(vector<Point_Type>::iterator it=critical_points.begin(); it!=critical_points.end(); ++it)
    {
        switch(*it)
        {
        case Point_Type::REGULAR:
            num_reg++;
            break;
        case Point_Type::MINIMUM:
            num_min++;
            break;
        case Point_Type::SADDLE:
            num_saddle++;
            break;
        case Point_Type::MULTIPLE_SADDLE:
            num_multisaddle++;
            break;
        case Point_Type::MAXIMUM:
            num_max++;
            break;
        }
    }
    cerr<<"[STAT] Extracted Critical points"<<endl;
    cerr<<"   mesh vertices: "<<critical_points.size()<<endl;
    cerr<<"   regular_points: "<<num_reg<<" -- minima: "<<num_min<<endl;
    cerr<<"   saddles: "<<num_saddle<<" -- multi-saddles: "<<num_multisaddle<<endl;
    cerr<<"   maxima: "<<num_max<<endl;
}
