/*
    This file is part of the Terrain Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)
               Paola Magillo (paola.magillo@unige.it)

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
    along with the triangle Trees library.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "c_curvature.h"
#include "utilities/timer.h"

void C_Curvature::curvature_leaf(Node_V &n, Mesh &mesh)
{
    /// if there are no vertices in the leaf we have nothing to do..
    if(!n.indexes_vertices())
        return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();
    itype v_range = v_end - v_start;

    leaf_VT vts(v_range,VT());
    leaf_VV vvs(v_range,VV());
    dvect v_aux(v_range,0.0);
    vector<Vertex> v_norm(v_range,Vertex());
    Vertex tnorm = Vertex();

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);

        Geometry_Curvature::triangle_normal(t,tnorm,mesh);

        for(int j=0; j<t.vertices_num(); j++)
        {
            itype real_index = t.TV(j);
            if(n.indexes_vertex(real_index))
            {
                this->update_local_structures(real_index,j,t,*t_id,tnorm,vts,vvs,v_aux,v_norm,v_start,/*is_v_border,*/mesh);
            }
        }
    }

    this->finalize_local_curvatures(v_start,vts,vvs,v_aux,v_norm,mesh);
}

void C_Curvature::curvature_leaf(Node_T &n, Box &n_dom, Mesh &mesh)
{
    itype v_start;
    itype v_end;

    n.get_v_range(v_start,v_end,n_dom,mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;

    leaf_VT vts(v_end-v_start,VT());
    leaf_VV vvs(v_end-v_start,VV());
    dvect v_aux(v_end-v_start,0.0);
    vector<Vertex> v_norm(v_end-v_start,Vertex());
    Vertex tnorm = Vertex();

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);

        Geometry_Curvature::triangle_normal(t,tnorm,mesh);

        for(int v_pos=0; v_pos<t.vertices_num(); v_pos++)
        {
            itype v_id = t.TV(v_pos);
            if(n.indexes_vertex(v_start,v_end,v_id))
            {
                this->update_local_structures(v_id,v_pos,t,*t_id,tnorm,vts,vvs,v_aux,v_norm,v_start,/*is_v_border,*/mesh);
            }
        }
    }

    this->finalize_local_curvatures(v_start,vts,vvs,v_aux,v_norm,mesh);
}

void C_Curvature::update_local_structures(itype v_id, itype v_pos, Triangle &t, itype t_id, Vertex &tnorm, leaf_VT &vts, leaf_VV &vvs, dvect &v_aux,
                                          vector<Vertex> &v_norm, itype v_start,/*boost::dynamic_bitset<> &is_v_border,*/ Mesh &mesh)
{
    int local_index = v_id-v_start;
    //VT
    vts[local_index].push_back(t_id);

    itype v1 = t.TV((v_pos+1)%t.vertices_num());
    itype v2 = t.TV((v_pos+2)%t.vertices_num());

    //VV
    vvs[local_index].insert(v1);
    vvs[local_index].insert(v2);

    coord_type ang = acos(Geometry_Curvature::cos_angle(mesh.get_vertex(v1),mesh.get_vertex(v_id),mesh.get_vertex(v2)));
    v_aux[local_index] += ang;

    for(int i=0; i <= tnorm.get_dimension(); i++)
        v_norm[local_index].set_c(i,v_norm[local_index].get_c(i) + (tnorm.get_c(i) * ang));
}

void C_Curvature::finalize_local_curvatures(itype v_start, leaf_VT &vts, leaf_VV &vvs, dvect &v_aux, vector<Vertex> &v_norm, Mesh &mesh)
{
    coord_type l, cc, min_cur, max_cur;

    /// differenza rispetto all'algo di Paola
    /// non visito i triangoli nell'intorno di un vertice interno, ma direttamente

    for(utype v=0; v<vvs.size(); v++)
    {
        l=0;
        min_cur=INFINITY;
        max_cur=-INFINITY;

        itype v_index = v + v_start;
        VT &vt = vts[v];
        VV &vv = vvs[v];

        if (vt.size()==(utype)1)
        {
            mesh.get_vertex(v_index).add_field(0.0); //one incident, we say that v is flat vertex
            continue;
        }

        for(int i=0; i <= v_norm[v].get_dimension(); i++)
        {
            v_norm[v].set_c(i,(v_norm[v].get_c(i) / (vt.size()*v_aux[v])));
            l += v_norm[v].get_c(i)*v_norm[v].get_c(i);
        }
        l = sqrt(l);

        for(int i=0; i <= v_norm[v].get_dimension(); i++)
            v_norm[v].set_c(i,(v_norm[v].get_c(i) / l));

        v_aux[v] = 0; // reset the v_aux variable

        for(auto vid : vv)
        {
            cc = compute_curve_and_curvature(v_index,mesh.get_vertex(v_index),vt,v_norm[v], mesh.get_vertex(vid),mesh);

            if (cc<min_cur)
                min_cur = cc;
            if (cc>max_cur)
                max_cur = cc;

            v_aux[v] += cc; // really used only if (take_mean_of_all)
        }

        if (type==MEAN)
        {
            if (take_mean_of_all)
                mesh.get_vertex(v_index).add_field(v_aux[v]/vv.size());
            else
                mesh.get_vertex(v_index).add_field((min_cur+max_cur)/2.0);
        }
        //gaussian Ccurvature
        else if (type==GAUSS)
            mesh.get_vertex(v_index).add_field(min_cur*max_cur);
    }
}

//compute curvature of polyline segment p1 v2 p3 and set its sign
//according to the direction of surface normal v2norm at vertex v2
coord_type C_Curvature::wedge_curvature(Vertex &p1, Vertex &v2, Vertex &p3, Vertex &v2norm)
{
    coord_type cosinus, curva;

    cosinus = Geometry_Curvature::cos_angle(p1,v2,p3);
    if (cosinus<=-1.0) curva = 0.0;
    else curva = PI - acos(cosinus);

    if (curva==0.0)
        return curva;

    Vertex v2_plus_norm = Vertex(v2.get_x()+v2norm.get_x(),v2.get_y()+v2norm.get_y(),v2.get_z()+v2norm.get_z());
    coord_type cos2, cos3, ang1, ang2, ang3;

    if (cosinus<=-1.0)
        ang1 = PI;
    else
        ang1 = acos(cosinus);
    cos2 = Geometry_Curvature::cos_angle(p1,v2,v2_plus_norm);
    if (cos2<=-1.0)
        ang2 = PI;
    else
        ang2 = acos(cos2);
    cos3 = Geometry_Curvature::cos_angle(p3,v2,v2_plus_norm);
    if (cos3<=-1.0)
        ang3 = PI;
    else
        ang3 = acos(cos3);
    // Note that ang1 is the angle <180 among the two angle formed by p1,v2,p3
    // so if surface normal lies inside this angle then surface is concave.
    // Therefore: change curvature sign iff ang1 == ang2+ang3,
    // but numerical errors may be present in the angles,
    // so do the following check:
    if ( ((ang2+ang3)<PI) && (ang2<ang1) && (ang3<ang1) )
        return -curva;
    else
        return curva;
}

//Given vertex v, its normal vnorm, and a point w, compute the curvature
//in v of the line obtained by intersecting the surface with the plane
//through v, vnorm, and w.
coord_type C_Curvature::compute_curve_and_curvature(itype v_id, Vertex &v, ivect &vt, Vertex &vnorm, Vertex &w, Mesh &mesh)
{
    dvect abc; abc.assign(3,0);// coefficients of plane through v,w and parallel to normal

    abc[0] =  +vnorm.get_z()*(v.get_y()-w.get_y()) -vnorm.get_y()*(v.get_z()-w.get_z());
    abc[1] =  -vnorm.get_z()*(v.get_x()-w.get_x()) +vnorm.get_x()*(v.get_z()-w.get_z());
    abc[2] =  -vnorm.get_x()*(v.get_y()-w.get_y()) +vnorm.get_y()*(v.get_x()-w.get_x());
    // equation plane:
    // a(x-vx)+b(y-vy)+c(z-vz)=0 dove v=(vx,vy,vz)

    Vertex w1 = Geometry_Curvature::find_plane_intersection(v_id, v, vt, vnorm, w, abc, mesh);

    //take curvature of wedge w v w1
    return wedge_curvature(w,v,w1, vnorm);
}

//constructor: read mesh and compute curvature for each vertex
C_Curvature :: C_Curvature(CCurvatureType tp, bool mean_of_all)
    : Abstract_Curvature()
{
    type = tp;
    take_mean_of_all = mean_of_all;
}
