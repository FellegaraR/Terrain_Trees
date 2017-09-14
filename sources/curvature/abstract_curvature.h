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

#ifndef ABSTRACT_CURVATURE
#define ABSTRACT_CURVATURE

#include "geometry/geometry_curvature.h"
#include "terrain_trees/node_v.h"
#include "terrain_trees/node_t.h"
#include "terrain_trees/spatial_subdivision.h"

class Abstract_Curvature
{
public:
    //compute for all vertices
    /**
     * @brief compute
     * @param n
     * @param mesh
     * @param division
     */
    void compute(Node_V &n, Mesh &mesh, Spatial_Subdivision &division);
    /**
     * @brief compute
     * @param n
     * @param domain
     * @param level
     * @param mesh
     * @param division
     */
    void compute(Node_T &n, Box &domain, int level, Mesh &mesh, Spatial_Subdivision &division);

    inline void print_curvature_stats(Mesh &mesh, int c_pos)
    {
        coord_type min=INFINITY, max=-INFINITY, avg=0;
        for(itype v=1; v<=mesh.get_vertices_num(); v++)
        {
            coord_type c = mesh.get_vertex(v).get_field(c_pos);

            if(c < min)
                min = c;
            if(c > max)
                max = c;
            avg += c;
        }
        cerr<<"[STATS] curvature min: "<<min<<" avg: "<<avg/(coord_type)mesh.get_vertices_num()<<" max: "<<max<<endl;
    }

protected:
    /**
     * @brief Abstract_Curvature
     * @param numVertex
     */
    Abstract_Curvature() { }
    /**
     * @brief curvature_leaf
     * @param n
     * @param mesh
     */
    virtual void curvature_leaf(Node_V& n, Mesh& mesh)=0;
    /**
     * @brief curvature_leaf
     * @param n
     * @param mesh
     */
    virtual void curvature_leaf(Node_T& n, Box &domain, Mesh& mesh)=0;
};

#endif // ABSTRACT_CURVATURE

