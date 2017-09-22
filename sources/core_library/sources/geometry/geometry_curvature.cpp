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

#include "geometry_curvature.h"

//compute cosinus of angle formed by 3 vertices
coord_type Geometry_Curvature::cos_angle(Vertex &v1, Vertex &v2, Vertex &v3)
{
    coord_type norm21=v2.norm(v1);
    coord_type norm23=v2.norm(v3);
    coord_type product=v2.scalar_product(v1,v3);
    coord_type costeta=product/(norm21*norm23);
    if (costeta>1.0) costeta = 1.0;
    if (costeta<-1.0) costeta = -1.0;
    return costeta;
}

//compute total area of triangles incident in v
coord_type Geometry_Curvature::fan_area(ivect& vt, Mesh& mesh)
{
    coord_type a = 0.0;

    for(unsigned int i=0;i<vt.size();i++)
    {
        a += triangle_area(mesh.get_triangle(vt[i]),mesh);
    }

    return a;
}

//compute triangle area
coord_type Geometry_Curvature::triangle_area(Triangle& t, Mesh &mesh)
{
    coord_type prodscaluv;
    coord_type normau;
    coord_type normav;
    coord_type cosalpha;
    coord_type senalpha;

    Vertex &v1 = mesh.get_vertex(t.TV(0));
    Vertex &v2 = mesh.get_vertex(t.TV(1));
    Vertex &v3 = mesh.get_vertex(t.TV(2));

    prodscaluv=v1.scalar_product(v2,v3);
    normau=v1.norm(v2);
    normav=v1.norm(v3);
    cosalpha = prodscaluv / (normau*normav);
    senalpha=sqrt(1-(cosalpha*cosalpha));
    if(std::isnan(senalpha)) senalpha=0.0001;

    return (normau*normav*senalpha)/2;
}

void Geometry_Curvature::triangle_normal(Triangle &t, Vertex &tnorm, Mesh &mesh)
{
    Vertex &v1 = mesh.get_vertex(t.TV(0));
    Vertex &v2 = mesh.get_vertex(t.TV(1));
    Vertex &v3 = mesh.get_vertex(t.TV(2));

    coord_type a[3], b[3]/*, x,y,z*/;
    coord_type norm=0;

    for(int i=0; i <= v1.get_dimension(); i++)
    {
        a[i] = v1.get_c(i)-v2.get_c(i);
        b[i] = v1.get_c(i)-v3.get_c(i);
    }

    tnorm.set_c(0,a[1]*b[2] - a[2]*b[1]);
    tnorm.set_c(1,a[2]*b[0] - a[0]*b[2]);
    tnorm.set_c(2,a[0]*b[1] - a[1]*b[0]);

    for(int i=0; i <= tnorm.get_dimension(); i++)
        norm += tnorm.get_c(i)*tnorm.get_c(i);
    norm = sqrt(norm);

    for(int i=0; i <= tnorm.get_dimension(); i++)
        tnorm.set_c(i,tnorm.get_c(i)/norm);
}

//test if segment v1-v2 intersect plane of equation
// a(x-vx)+b(y-vy)+c(z-vz)=0 where v=(vx,vy,vz), three possible results:
// 0 : no intersection
// 1 : intersection is v1
// 2 : intersection is v2
// 3 : proper intersection
int Geometry_Curvature::intersect_plane(Vertex &v1, Vertex &v2, dvect &abc, Vertex &v)
{
    coord_type res1=0, res2=0;
    for(int i=0; i <= v1.get_dimension(); i++)
    {
        res1 += abc[i]*(v1.get_c(i)-v.get_c(i));
        res2 += abc[i]*(v2.get_c(i)-v.get_c(i));
    }

    if ((-TOLER<=res1)&&(res1<=TOLER)) return 1;
    if ((-TOLER<=res2)&&(res2<=TOLER)) return 2;
    if ( (res1>0.0)&&(res2<0.0) ) return 3;
    if ( (res2>0.0)&&(res1<0.0) ) return 3;
    return 0;
}

int Geometry_Curvature::same_point(Vertex &p1, Vertex &p2)
{
    for(int i=0; i <= p1.get_dimension(); i++) /// NOTA <=
        if (fabs(p1.get_c(i)-p2.get_c(i))>SMALL_TOLER)
            return 0;
    return 1;
}

//Find the triangle, around vertex v, that is intersected by the plane through
//v, the normal vnorm to v, and point w, the equation of such plane is
//a(x-vx)+b(y-vy)+c(z-vz)=0 dove v=(vx,vy,vz):
//compute the intersection point and return it.
Vertex Geometry_Curvature::find_plane_intersection(itype v_id, Vertex &v, ivect &vt, Vertex &vnorm, Vertex &w,
                                                   dvect &abc, Mesh &mesh)
{
    bool found = false;
    itype v1_id, v2_id;
    Vertex w1;

    for(ivect_iter it=vt.begin(); it!=vt.end(); ++it)
    {
        Triangle &t = mesh.get_triangle(*it);
        int v_pos = t.vertex_index(v_id);

        v1_id = t.TV((v_pos+1)%t.vertices_num());
        v2_id = t.TV((v_pos+2)%t.vertices_num());
        Vertex &v1 = mesh.get_vertex(v1_id);
        Vertex &v2 = mesh.get_vertex(v2_id);

        if ( !Geometry_Curvature::same_point(v1,w) && !Geometry_Curvature::same_point(v2,w) )
        {
            switch(Geometry_Curvature::intersect_plane(v1,v2,abc,v))
            {
            case 1: // intersection is v1
                w1 = v1;
                found = true;
                break;
            case 2: // intersection is v2
                w1 = v2;
                found = true;
                break;
            case 3: // proper intersection, compute intersection                
                coord_type d=0,n=0;
                for(int i=0; i <= v.get_dimension(); i++)
                {
                    d += -abc[i]*v1.get_c(i) + abc[i]*v.get_c(i);
                    n +=  abc[i]*v2.get_c(i) - abc[i]*v1.get_c(i);
                }
                coord_type s = d / n;
                for(int i=0; i <= v1.get_dimension(); i++) /// NOTA <=
                    w1.set_c(i,s*(v2.get_c(i)-v1.get_c(i)) + v1.get_c(i));
                found = true;
                break;
            }

            // se ho trovato w stesso devo cercare l'altro punto di intersezione
            if (found && same_point(w1,w))
                found = false;
        }

        if(found)
            break;
    }

    if(!found)
        // v is a boundary vertex, we have visited all triangles
        // and not found an edge which intersects...
    {
        // simulate a vertical wall
        for(int i=0; i <= v.get_dimension(); i++) /// NOTA <=
            w1.set_c(i,v.get_c(i)-vnorm.get_c(i));
    }

    return w1;
}
