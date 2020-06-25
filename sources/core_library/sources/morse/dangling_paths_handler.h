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

#ifndef DANGLING_PATHS_HANDLER_H
#define DANGLING_PATHS_HANDLER_H

#include "forman_gradient.h"
#include "terrain_trees/node_v.h"
#include "terrain_trees/node_t.h"
#include "terrain_trees/spatial_subdivision.h"

class Dangling_Paths_Handler
{
public:
    static Node_V& find_leaf(itype v, Node_V &n, Spatial_Subdivision &division, itype &key);
    static Node_T& find_leaf(Vertex &v, Node_T &n, Box &n_dom, int n_level, Mesh &mesh, Spatial_Subdivision &division, itype &other_v_start, itype &other_v_end);

    static void add_dangling_path(Node_V &root, Spatial_Subdivision &division, leaves_2_desc_map &map, itype t_id, Triangle& t, itype label);
    static void add_dangling_path(Node_V &root, Spatial_Subdivision &division, leaves_1_desc_map &map, Edge *edge, itype label);
    static void add_dangling_path(Node_V &root, Spatial_Subdivision &division, leaves_2_asc_map &map, itype v_id, itype current_v, itype label);
    static void add_dangling_path(Node_V &root, Spatial_Subdivision &division, leaves_asc_map &map,  itype max_v_id, itype t_id, short  e_pos, itype label);
    static void add_dangling_path(Node_V &root, Spatial_Subdivision &division, leaves_asc_map &map, itype max_v_id, pair<itype,short> &tf, itype label);

    ///MIG - ascending_1
    static void add_dangling_path(Node_V &root, Spatial_Subdivision &division, leaves_1_asc_mig_map &map,
                           itype max_v_id, pair<itype,short> &tf, itype label, iNode *saddle_node);
    static void add_dangling_path(Node_T &root, Mesh &mesh, Spatial_Subdivision &division, leaves_1_asc_mig_map &map,
                           itype max_v_id, pair<itype,short> &tf, itype label, iNode *saddle_node);
    ///MIG - descending_1
    static void add_dangling_path(Node_V &root, Spatial_Subdivision &division, leaves_1_desc_mig_map &map, Edge *edge, itype last_v, short label, iNode *saddle_node);
    static void add_dangling_path(Node_T &root, Mesh &mesh, Spatial_Subdivision &division, leaves_1_desc_mig_map &map, Edge *edge, itype last_v, short label, iNode *saddle_node);

private:
    Dangling_Paths_Handler() {}
};

#endif // DANGLING_PATHS_HANDLER_H
