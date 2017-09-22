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

#ifndef TREE_H
#define	TREE_H

#include "basic_types/mesh.h"
#include "basic_types/soup.h"
#include "spatial_subdivision.h"

///A super-class not instantiable representing a generic tree
template<class N> class Tree
{
public:
    ///A public method that returns the mesh associated to the tree
    /*!
     * \return a Mesh& variable, representing the mesh
     */
    inline Mesh& get_mesh() { return this->mesh; }
    ///A public method that returns the root node of the tree
    /*!
     * \return a N& variable, representing the root
     */
    inline N& get_root() { return this->root; }
    ///A public method that returns the spatial decomposition of the tree
    /*!
     * \return a Spatial_Decomposition& variable, representing the spatial decomposition of the tree
     */
    inline Spatial_Subdivision& get_subdivision() { return this->subdivision; }
    ///A public pure virtual method, implemented by the heirs class, that builds the tree
    virtual void build_tree()=0;    
    /*!
     * \brief A public pure virtual method, implemented by the heirs class, that builds the tree from a soup of triangles
     * \param soup
     */
    virtual void build_tree(Soup &soup)=0;
    ///A public pure virtual method, implemented by the heirs class, that builds the tree from a points cloud
    virtual void build_tree_from_cloud(vertex_multifield &multifield)=0;


    /**
     * @brief A public method that initialize the array of leaf blocks references in order to enable the parallel visit of leaf blocks
     *
     * @param n a N& argument representing the current node
     *
     * NOTA: this procedure must be called after loading/building a Stellar tree
     */
    void init_leaves_list(N& n);
    /**
     * @brief A public method that returns the leaf block at i-th position in the leaves array
     *
     * @param i an integer argument representing the position index of the leaf
     * @return N* a pointer reference to the leaf block
     */
    inline N* get_leaf(unsigned i) { return this->leaves[i]; }
    /**
     * @brief A public method that returns the number of leaf blocks in the tree
     *
     * @return unsigned
     */
    inline unsigned get_leaves_number() { return this->leaves.size(); }
    
protected:
    ///A constructor method
    Tree() { }
    ///A constructor method
    Tree(int sons_num)
    {
        this->subdivision = Spatial_Subdivision(sons_num);
    }
    ///A constructor method
    Tree(const Tree& orig)
    {
        this->subdivision = orig.subdivision;
        this->mesh = orig.mesh;
        this->root = orig.root;
    }
    ///A destructor method
    virtual ~Tree() {}

    ///A protected variable representing the mesh associated to the tree
    Mesh mesh;
    ///A protected variable representing the root node of the tree
    N root;
    ///A protected variable representing the division type of the tree
    Spatial_Subdivision subdivision;
    ///A protected variable that contains the list of pointers to the leaf block of the tree. This enables the PARALLELIZATION of the index lookup
    vector<N*> leaves;

    ///A protected pure virtual method, implemented by the heirs class, that split a node, creating the sons node, following the current division type
    /*!
     * \param n a N& argument, represents the node to split
     * \param domain a Box& argument, represents the node domain
     * \param level an integer argument, representing the node level in the hierarchy
     */
    virtual void split(N& n, Box& domain, int level)=0;

private:

};

template<class N> void Tree<N>::init_leaves_list(N& n)
{
    if(!n.is_leaf())
    {
        for(typename N::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
            {
                if((*it)->is_leaf())
                {
                    this->leaves.push_back(*it);
                }
                else
                {
                    this->init_leaves_list(**it);
                }
            }
        }
    }
}

#endif	/* TREE_H */

