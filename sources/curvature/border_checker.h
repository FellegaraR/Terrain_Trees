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

#ifndef BORDERCHECKER_H
#define BORDERCHECKER_H

#include <iostream>
#include <algorithm>
#include <map>

#include "utilities/sorting.h"
#include "terrain_trees/node_v.h"
#include "terrain_trees/node_t.h"
#include "terrain_trees/spatial_subdivision.h"

/**
 * @brief The Border_Checker class exploits the borders of a triangle mesh indexed by a Terrain tree
 * The mesh borders are needed when computing the Concentrated Curvature.
 */
class Border_Checker
{
public:
    Border_Checker() {}
    /**
     * @brief A public procedure that visits recursively a tree and exploit the mesh borders.
     * This version is compatible with PMR-T trees on which the spatial coherence has been exploited.
     *
     * @param n a Node_T& parameter representing the current node
     * @param dom a Box& parameter representing the current domain of node n
     * @param level an integer argument representing the level of n in the hierarchy
     * @param mesh a Mesh& parameter representing the triangle mesh
     * @param division a Spatial_Subdivision& parameter representing the spatial subdivision of the tree
     */
    void compute_borders(Node_T &n, Box &dom, int level, Mesh &mesh, Spatial_Subdivision &division);
    /**
     * @brief A public procedure that visits recursively a tree and exploit the mesh borders.
     * This version requires a tree on which the spatial coherence has been exploited, and is compatible with PR-T and PM-T trees.
     *
     * @param n a Node_V& parameter representing the current node
     * @param dom a Box& parameter representing the current domain of node n
     * @param level an integer argument representing the level of n in the hierarchy
     * @param mesh a Mesh& parameter representing the triangle mesh
     * @param division a Spatial_Subdivision& parameter representing the spatial subdivision of the tree
     */
    void compute_borders(Node_V &n, Box &dom, int level, Mesh &mesh, Spatial_Subdivision &division);

private:
    /**
     * @brief A private procedure that exploits the mesh border in a leaf block.
     * This version requires a tree on which the spatial coherence has been exploited, and is compatible with P-Ttrees and PT-Ttrees.
     *
     * @param n a Node_V& parameter representing the current node
     * @param mesh a Mesh& parameter representing the triangle mesh
     */
    void compute_borders(Node_V& n, Mesh& mesh);
    /**
     * @brief A private procedure that exploits the mesh border in a leaf block.
     * This version is compatible with T-Ttrees and RT-Ttrees on which the spatial coherence has been exploited
     * as well as on Terrain Trees on which the spatial coherence is not exploited.
     *
     * @param n a Node_T& parameter representing the current node
     * @param dom a Box& parameter representing the current domain of node n
     * @param mesh a Mesh& parameter representing the triangle mesh
     */
    void compute_borders(Node_T& n, Box &n_dom, Mesh& mesh);
    /**
     * @brief A private procedure that sets the border of the sub-mesh indexed by the leaf block
     *
     * @param edges an array containing the tuples composed by the 2 indices forming and and the triangle in its co-boundary
     * @param mesh a Mesh& parameter representing the triangle mesh
     * @return true if we change the borders, false otherwise
     */
    bool set_borders(vector<edge_triangle_tuple> &edges, Mesh& mesh);
    /**
     * @brief A private procedure that extracts the two edges incident in a vertex in the boundary of a triangle
     *
     * @param t a Triangle& parameter, representing the current triangle
     * @param t_id an itype parameter, representing the position index of t
     * @param v_pos an integer parameter representing the vertex position index in the boundary of t
     * @param edges an array containing the edges incident in the vertex
     */
    void get_incident_edges(Triangle &t, itype t_id, int v_pos, vector<edge_triangle_tuple> &edges);
};



#endif // BORDERCHECKER_H
