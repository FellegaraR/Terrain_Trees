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
    along with the triangle Trees library.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef GEOMETRY_CURVATURE_H
#define GEOMETRY_CURVATURE_H

#include "geometry.h"
#include "basic_types/basic_wrappers.h"
#include "basic_types/mesh.h"

#define TOLER (0.3e-6)
#define SMALL_TOLER (0.3e-10)

class Geometry_Curvature
{
public:
    Geometry_Curvature() {}

    //compute cosinus of angle formed by 3 vertices
    static coord_type cos_angle(Vertex& v1, Vertex& v2, Vertex& v3);
    //compute total area of triangles incident in v
    static coord_type fan_area(ivect &vt, Mesh &mesh);
    //compute triangle area
    static coord_type triangle_area(Triangle& t, Mesh& mesh);
    //compute triangle normal
    static void triangle_normal(Triangle &t, Vertex &tnorm, Mesh &mesh);
    //test if segment v1-v2 intersect plane of equation:
    // a(x-vx)+b(y-vy)+c(z-vz)=0 where v=(vx,vy,vz)
    static int intersect_plane(Vertex &v1, Vertex &v2, dvect &abc, Vertex &v);
    //Find the triangle, around vertex v, that is intersected by the plane through
    //v, the normal vnorm to v, and point w, the equation of such plane is
    //a(x-vx)+b(y-vy)+c(z-vz)=0 dove v=(vx,vy,vz):
    //compute the intersection point and return it.
    static Vertex find_plane_intersection(itype v_id, Vertex &v, ivect &vt, Vertex &vnorm, Vertex &w,
                                          dvect &abc, Mesh &mesh);
    static int same_point(Vertex &p1, Vertex &p2);
};

#endif // GEOMETRY_CURVATURE_H
