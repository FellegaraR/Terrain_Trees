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

#include "sorting.h"
#include <algorithm>

void sorting_vertices(vector<vertex_triangle_pair> &vert_vec, const ivect &triangles, Mesh &mesh)
{
    int k=0;
    //create the array of pairs vertex-triangle
    for(ivect_citer t_id=triangles.begin(); t_id!=triangles.end(); ++t_id)
    {
        Triangle &t = mesh.get_triangle(*t_id);
        for(int i=0;i<t.vertices_num();i++)
        {
            vert_vec[k] = vertex_triangle_pair(t.TV(i),*t_id);
            k++;
        }
    }
    //order the array
    std::sort(vert_vec.begin(),vert_vec.end());
}

void sort_edges_tuples(vector< edge_triangle_tuple >& edges) { std::sort(edges.begin(),edges.end()); }
