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

#ifndef GEOMETRY_SLOPE_H
#define GEOMETRY_SLOPE_H

//#include "geometry.h"
#include "basic_types/basic_wrappers.h"
#include "basic_types/mesh.h"

class Geometry_Slope
{
public:
    static coord_type compute_triangle_slope(Triangle &t, Mesh &mesh);
    
    static coord_type compute_triangle_slope_song(Triangle &t, Mesh &mesh);

    static coord_type compute_edge_slope(ivect &e, Mesh &mesh);

private:
    Geometry_Slope() {}
};

#endif // GEOMETRY_SLOPE_H
