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

#ifndef NODE_T_H
#define NODE_T_H

#include "terrain_trees/node.h"

/**
 * @brief The Node_T class extends the Node class, providing an instantiable definition of it.
 * This class is used by the PMR-T trees
 */
class Node_T : public Node<Node_T>
{
public:
    ///A constructor
    Node_T() : Node<Node_T>() { }
    ///A copy-constructor
    Node_T(const Node_T& orig) : Node<Node_T>(orig) { }
    ///A destructor
    virtual ~Node_T() {}
    /**
     * @brief operator <<
     * @param out
     * @param p
     * @return
     */
    friend std::ostream& operator<<(std::ostream& out, const Node_T& p)
    {
        if(p.is_leaf())
            out <<"Leaf[t->"<< p.get_real_t_array_size() <<"]";
        else
            out <<"Node";
        return out;
    }

    /**
     * @brief A public method that returns the range of vertices completely contained into the leaf
     * NOTA: this method only works after the reordering of vertices array
     *
     * @param v_start an integer that it will contains the first vertex indexed by the leaf
     * @param v_end an integer that it will contains the first vertex outside the leaf
     * @param dom a Box& representing the node domain
     * @param mesh a Mesh& representing the triangle mesh
     */
    void get_v_range(itype &v_start, itype &v_end, Box& dom, Mesh& mesh);
    /**
     * @brief A public method that checks if a vertex is indexed by the current node
     * NOTA: this method is a wrapper thought for T-Ttrees and RT-Ttrees that do not explicitly encode the vertices
     * within each node. For P-Ttrees and PT-Ttrees a wrapper function has been defined
     *
     * @param v_start an integer that it will contains the first vertex indexed by the leaf
     * @param v_end an integer that it will contains the first vertex outside the leaf
     * @param v_id an integer representing the position index of a vertex
     * @return true if v_id is into the leaf, false otherwise
     */
    inline bool indexes_vertex(itype v_start, itype v_end, itype v_id) { return (v_id >= v_start && v_id < v_end); }

    /**
     * @brief A public method that extracts the Vertex-Triangle (VT) relations for the vertices indexed in the current block
     * @param all_vt a leaf_VT variable, that encodes the VT relations
     * @param dom a Box& representing the node domain
     * @param mesh a Mesh& variable representing the triangle mesh
     */
    void get_VT(leaf_VT &all_vt, Box &dom, Mesh &mesh);
    /**
     * @brief A public method that extracts the Vertex-Triangle (VT) relations for the vertices indexed in the current block
     * @param all_vt a leaf_VT variable, that encodes the VT relations
     * @param v_start an integer that it will contains the first vertex indexed by the leaf
     * @param v_end an integer that it will contains the first vertex outside the leaf
     * @param mesh a Mesh& variable representing the triangle mesh
     */
    void get_VT(leaf_VT &all_vt, itype v_start, itype v_end, Mesh &mesh);
    /**
     * @brief A public method that extracts the Vertex-Vertex (VV) relations for the vertices indexed in the current block
     * @param all_vv a leaf_VV variable, that encodes the VV relations
     * @param dom a Box& representing the node domain
     * @param mesh a Mesh& variable representing the triangle mesh
     */
    void get_VV(leaf_VV &all_vv, Box &dom, Mesh& mesh);
    /**
     * @brief A public method that extracts the Vertex-Vertex (VV) relations for the vertices indexed in the current block
     * @param all_vv a leaf_VV variable, that encodes the VV relations
     * @param v_start an integer that it will contains the first vertex indexed by the leaf
     * @param v_end an integer that it will contains the first vertex outside the leaf
     * @param mesh a Mesh& variable representing the triangle mesh
     */
    void get_VV(leaf_VV &all_vv, itype v_start, itype v_end, Mesh& mesh);
    void get_VV_vector(leaf_VV_vec &all_vv, Box &dom, Mesh& mesh);
    void get_VV_vector(leaf_VV_vec &all_vv, itype v_start, itype v_end, Mesh& mesh);
    /**
     * @brief A public method that extracts the Vertex-Triangle (VT) and Vertex-Vertex (VV) relations for the vertices indexed in the current block
     * @param all_vv a leaf_VV variable, that encodes the VV relations
     * @param all_vt a leaf_VT variable, that encodes the VT relations
     * @param v_start an integer that it will contains the first vertex indexed by the leaf
     * @param v_end an integer that it will contains the first vertex outside the leaf
     * @param mesh a Mesh& variable representing the triangle mesh
     */
    void get_VV_VT(leaf_VV &all_vv, leaf_VT &all_vt, itype v_start, itype v_end, Mesh& mesh);
    /**
     * @brief A public method that extracts the Edge-Triangle (ET) relations for the edges indexed in the current block
     * Notice that an edge is considered processed only if the current leaf block indexes the extreme with the higher position index
     * in this way each edge is processed once during each traversal of the tree
     *
     * @param ets a leaf_ET variable, that encodes the ET relations
     * @param v_start an integer that it will contains the first vertex indexed by the leaf
     * @param v_end an integer that it will contains the first vertex outside the leaf
     * @param mesh a Mesh& variable representing the triangle mesh
     */
    void get_ET(leaf_ET &ets, itype v_start, itype v_end, Mesh &mesh);

protected:
};

#endif // NODE_T_H
