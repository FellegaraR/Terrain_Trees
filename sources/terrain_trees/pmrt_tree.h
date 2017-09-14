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

#ifndef RT_TREE_H
#define	RT_TREE_H

#include "tree.h"
#include "node_t.h"
#include <geometry/geometry_wrapper.h>

///An inner-class, implementing Tree, that represents a tree which uses the PMRT-T criterion
class PMRT_Tree : public Tree<Node_T>
{
public:
    ///A constructor method
    PMRT_Tree() { }
    ///A constructor method
    /*!
     * \param triangles_per_leaf an integer argument, represents the maximum number of triangles admitted for a tree node
     */
    PMRT_Tree(int triangles_per_leaf, int sons_num);
    ///A constructor method
    PMRT_Tree(const PMRT_Tree& orig) : Tree<Node_T>(orig) { this->triangles_threshold = orig.triangles_threshold; }
    ///A destructor method
    virtual ~PMRT_Tree() {}
    ///A public method that sets the maximum number of triangles admitted for a node
    /*!
     * \param maxT an integer argument, represents the maximum number of triangles
     */
    inline void set_triangles_threshold(int maxT) { this->triangles_threshold = maxT; }
    ///A public method that returns the maximum number of triangles admitted for a node
    /*!
     * \return an integer value representing the maximum number of triangles admitted
     */
    inline int get_triangles_threshold() { return this->triangles_threshold; }
    ///A public method that builds the tree
    /**
     * Only the triangles are inserted into the hierarchy.
     */
    void build_tree();
    /*!
     * \brief A public method that builds the tree from a soup of triangles
     * \param soup
     */
    void build_tree(Soup &)
    {
        std::cerr<<"[ERROR] build_tree: the tree generation from a soup of triangles is not supported on PMR-T Trees."<<endl;
    }
    ///A public method that builds the tree from a points cloud
    void build_tree_from_cloud(vertex_multifield &)
    {
        std::cerr<<"[ERROR] build_tree: the tree generation from a points cloud is not supported on PMR-T Trees."<<endl;
    }

private:
    ///A private variable representing the maximum number of triangles admitted for a node
    int triangles_threshold;
    ///A private method that adds a triangle to the tree structure
    /*!
     * This method checks if a triangle has a proper intersection with the node, and then insert the triangle into
     * the node list if the node is a leaf, otherwise if the node is internal a recursive call is activated.
     * If the maximum number of triangles for a node is exceeded, a splitting operation is activated only once.
     *
     * \param n a Node_T& argument, represents the node in which we try to insert the triangle
     * \param domain a Box& argument, represents the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * \param t an interger argument, represents the triangle to insert
     */
    void add_triangle(Node_T& n, Box& domain, int level, itype t);
    ///A private method that reinsert a triangle to the tree structure only once
    /*!
     * This method checks if a triangle has a proper intersection with the node, and then insert the triangle into
     * the node list if the node is a leaf. The block threshold is not checked in this method, thus, the every time
     * a split operation is called only on the splitted parent node.
     *
     * \param n a Node_T& argument, represents the node in which we try to insert the triangle
     * \param domain a Box& argument, represents the node domain
     * \param t an integer argument, represents the triangle to insert
     */
    void reinsert_triangle_once(Node_T& n, Box& domain, itype t);
    ///A protected method that split a node, creating the sons node, following the current division type
    /*!
     * This method also reinsert all the triangle, which are in the splitted node, into the sons node.
     * But unlike the other criteria it calls the addtriangle_once method, which doesn't cause a recursive insertion into the sons node.
     *
     * \param n a Node_T& argument, represents the node to split
     * \param domain a Box& argument, represents the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     */
    void split(Node_T& n, Box& domain, int level);

    ///A public method that checks if a node contains the maximum number of triangles admitted
    /*!
     * \param n a Node_T& argument, represents the node to check
     * \return a boolean value, true if the limit is exceeded, false otherwise
     */
    inline bool is_full(Node_T &n) { return (n.get_t_array_size() > this->triangles_threshold); }
};

#endif	/* RT_TREE_H */

