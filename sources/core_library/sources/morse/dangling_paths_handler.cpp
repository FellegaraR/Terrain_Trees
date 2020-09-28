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

#include "dangling_paths_handler.h"

Node_V& Dangling_Paths_Handler::find_leaf(itype v, Node_V &n, Spatial_Subdivision &division, itype &key)
{
    if (n.is_leaf())
    {
        key = (n.get_v_start()+n.get_v_end());
        return n;
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            //we visit a son only if it contains the vertex
            if(n.get_son(i)!=NULL && n.get_son(i)->indexes_vertex(v))
            {
                return Dangling_Paths_Handler::find_leaf(v,*n.get_son(i),division,key);
            }
        }
    }
}

Node_T& Dangling_Paths_Handler::find_leaf(Vertex &v, Node_T &n, Box &n_dom, int n_level, Mesh &mesh, Spatial_Subdivision &division,
                                          itype &other_v_start, itype &other_v_end)
{
    if (n.is_leaf())
    {
        n.get_v_range(other_v_start,other_v_end,n_dom,mesh); // we need to gather the vertices range..
        return n;
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(n_dom,n_level,i);
            int son_level = n_level +1;
            //we visit a son only if it contains the vertex
            if(n.get_son(i)!=NULL && son_dom.contains(v,mesh.get_domain().get_max()))
            {
                return Dangling_Paths_Handler::find_leaf(v,*n.get_son(i),son_dom,son_level,mesh,division,other_v_start,other_v_end);
            }
        }
    }
}

void Dangling_Paths_Handler::add_dangling_path(Node_V &root, Spatial_Subdivision &division, leaves_2_desc_map &map, itype t_id, Triangle &t, itype label)
{
    itype key = -1;
    itype max_index = t.maxindex();
    Dangling_Paths_Handler::find_leaf(max_index,root,division,key);

    leaves_2_desc_map::iterator iter = map.find(key);
    if(iter != map.end())
    {
        map_2paths &inside = iter->second;
        inside.insert((make_pair(t_id,label)));
    }
    else
    {
        map_2paths new_one;
        new_one.insert((make_pair(t_id,label)));
        map[key] = new_one;
    }
}

void Dangling_Paths_Handler::add_dangling_path(Node_V &root, Spatial_Subdivision &division, leaves_1_desc_map &map, Edge *edge, itype label)
{
    itype key = -1;
    Dangling_Paths_Handler::find_leaf(edge->maxindex(),root,division,key);
    leaves_1_desc_map::iterator iter = map.find(key);
    pair<itype,itype> coppia = make_pair(edge->EV(0),edge->EV(1));
    if(iter != map.end())
    {
        map_paths_pair &inside = iter->second;
        inside.insert(make_pair(coppia,label));
    }
    else
    {
        map_paths_pair new_one;
        new_one.insert((make_pair(coppia,label)));
        map[key] = new_one;
    }
}

void Dangling_Paths_Handler::add_dangling_path(Node_V &root, Spatial_Subdivision &division, leaves_2_asc_map &map, itype v_id, itype current_v, itype label)
{
    itype key = -1;
    Dangling_Paths_Handler::find_leaf(v_id,root,division,key);

    pair<itype,itype> coppia = make_pair(v_id,current_v);

    leaves_2_asc_map::iterator iter = map.find(key);
    if(iter != map.end())
    {
        map_paths_pair &inside = iter->second;
        inside.insert(make_pair(coppia,label));
    }
    else
    {
        map_paths_pair new_one;
        new_one.insert((make_pair(coppia,label)));
        map[key] = new_one;
    }
}

void Dangling_Paths_Handler::add_dangling_path(Node_V &root, Spatial_Subdivision &division, leaves_asc_map &map, itype max_v_id, itype t_id, short e_pos, itype label)
{
    itype key = -1;
    Dangling_Paths_Handler::find_leaf(max_v_id,root,division,key);
    pair<itype,short> coppia = make_pair(t_id,e_pos);
    triple new_triple = triple(coppia,label);

    leaves_asc_map::iterator iter = map.find(key);
    if(iter != map.end())
    {
        triple_set &inside = iter->second;
        inside.insert(new_triple);
    }
    else
    {
        triple_set new_one;
        new_one.insert(new_triple);
        map[key] = new_one;
    }
}

void Dangling_Paths_Handler::add_dangling_path(Node_V &root, Spatial_Subdivision &division, leaves_asc_map &map, itype max_v_id, pair<itype,short> &tf, itype label)
{
    itype key = -1;
    Dangling_Paths_Handler::find_leaf(max_v_id,root,division,key);
    triple new_triple = triple(tf,label);

    leaves_asc_map::iterator iter = map.find(key);
    if(iter != map.end())
    {
        triple_set &inside = iter->second;
        inside.insert(new_triple);
    }
    else
    {
        triple_set new_one;
        new_one.insert(new_triple);
        map[key] = new_one;
    }
}

void Dangling_Paths_Handler::add_dangling_path(Node_V &root, Spatial_Subdivision &division, leaves_1_asc_mig_map &map,
                                             itype max_v_id, pair<itype,short> &tf, itype label, iNode *saddle_node)
{
    itype key = -1;
    Dangling_Paths_Handler::find_leaf(max_v_id,root,division,key);
    quadruple new_quadruple = quadruple(tf,label,saddle_node);

    leaves_1_asc_mig_map::iterator iter = map.find(key);
    if(iter != map.end())
    {
        set_asc_quadruple &inside = iter->second;
        if(!inside.insert(new_quadruple).second)
        {
            cout<<"already into the dangling paths"<<endl;
            int a; cin>>a;
        }
    }
    else
    {
        set_asc_quadruple new_one;
        if(!new_one.insert(new_quadruple).second)
        {
            cout<<"already into the dangling paths"<<endl;
            int a; cin>>a;
        }
        map[key] = new_one;
    }
}

void Dangling_Paths_Handler::add_dangling_path(Node_T &root, Mesh &mesh, Spatial_Subdivision &division, leaves_1_asc_mig_map &map,
                                             itype max_v_id, pair<itype,short> &tf, itype label, iNode *saddle_node)
{
    itype other_v_start, other_v_end;
    Dangling_Paths_Handler::find_leaf(mesh.get_vertex(max_v_id),root,mesh.get_domain(),0,mesh,division,other_v_start,other_v_end);
    itype key = other_v_start+other_v_end;
    quadruple new_quadruple = quadruple(tf,label,saddle_node);

    leaves_1_asc_mig_map::iterator iter = map.find(key);
    if(iter != map.end())
    {
        set_asc_quadruple &inside = iter->second;
        if(!inside.insert(new_quadruple).second)
        {
            cout<<"already into the dangling paths"<<endl;
            int a; cin>>a;
        }
    }
    else
    {
        set_asc_quadruple new_one;
        if(!new_one.insert(new_quadruple).second)
        {
            cout<<"already into the dangling paths"<<endl;
            int a; cin>>a;
        }
        map[key] = new_one;
    }
}

void Dangling_Paths_Handler::add_dangling_path(Node_V &root, Spatial_Subdivision &division, leaves_1_desc_mig_map &map, Edge *edge,
                                               itype last_v, short label, iNode *saddle_node)
{
    itype key = -1;
    Dangling_Paths_Handler::find_leaf(edge->maxindex(),root,division,key);
    desc1_mig_quadruple q = desc1_mig_quadruple(edge,last_v,label,saddle_node);
    leaves_1_desc_mig_map::iterator iter = map.find(key);
    if(iter != map.end())
    {
        set_1paths &inside = iter->second;
        if(!inside.insert(q).second)
        {
            cout<<"already into the dangling paths"<<endl;
            int a; cin>>a;
        }
    }
    else
    {
        set_1paths new_one;
        if(!new_one.insert(q).second)
        {
            cout<<"already into the dangling paths"<<endl;
            int a; cin>>a;
        }
        map[key] = new_one;
    }
}

void Dangling_Paths_Handler::add_dangling_path(Node_T &root, Mesh &mesh, Spatial_Subdivision &division, leaves_1_desc_mig_map &map, Edge *edge,
                                               itype last_v, short label, iNode *saddle_node)
{
    itype other_v_start, other_v_end;
    Dangling_Paths_Handler::find_leaf(mesh.get_vertex(edge->maxindex()),root,mesh.get_domain(),0,mesh,division,other_v_start,other_v_end);
    itype key = other_v_start+other_v_end;
    desc1_mig_quadruple q = desc1_mig_quadruple(edge,last_v,label,saddle_node);
    leaves_1_desc_mig_map::iterator iter = map.find(key);
    if(iter != map.end())
    {
        set_1paths &inside = iter->second;
        if(!inside.insert(q).second)
        {
            cout<<"already into the dangling paths"<<endl;
            int a; cin>>a;
        }
    }
    else
    {
        set_1paths new_one;
        if(!new_one.insert(q).second)
        {
            cout<<"already into the dangling paths"<<endl;
            int a; cin>>a;
        }
        map[key] = new_one;
    }
}
