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

#ifndef _SPATIAL_DECOMPOSITION_H
#define	_SPATIAL_DECOMPOSITION_H

#include "basic_types/box.h"

///A super-class, not instantiable, representing a generic spatial decomposition of a tree
class Spatial_Subdivision
{
public:
    ///A constructor method
    Spatial_Subdivision() { sons = -1; }
    ///A constructor method
    /**
     * @brief Spatial_Decomposition
     * @param sn
     */
    Spatial_Subdivision(int sn) { sons = sn; }
    ///A copy-constructor method
    /**
     * @brief Spatial_Decomposition
     */
    Spatial_Subdivision(const Spatial_Subdivision &orig) { sons = orig.sons; }
    ///A destructor method
    virtual ~Spatial_Subdivision() {}
    ///Public method that returns the number of son nodes
    /*!
     * \return an integer value representing the number of sons
     */
    inline int son_number() { return sons; }
    ///Public method that computes the box domain of a node
    /*!
     * Depending from the
     *
     * \param parent_dom a Box& argument, representing the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * \param child_ind a integer argument, representing the index position of the son in the sub-tree
     * \return a Box value, representing the domain of the son
     */
    Box compute_domain(Box& parent_dom, int level, int child_ind);

private:
    int sons;

    ///Public method that computes the box domain of a son node following the quadtree decomposition
    /*!
     * \param parent_dom a Box& argument, representing the node domain
     * \param child_ind a integer argument, representing the son index position in the sub-tree
     * \return a Box value, representing the domain of the son node
     */
    Box compute_quad_domain(Box& parent_dom, int child_ind);
    ///Public method that computes the box domain of a son node following the kD-tree decomposition
    /*!
     * In this case the domain is computed considering the parent node level.
     * The subdivision splits alternately the two coordinate axes.
     *
     * \param parent_dom a Box& argument, representing the node domain
     * \param level an integer argument representing the level of n in the hierarchy
     * \param child_ind a integer argument, representing the son index position in the sub-tree
     * \return a Box value, representing the domain of the son node
     */
    Box compute_kd_domain(Box& parent_dom, int level, int child_ind);

};

#endif	/* _SPATIAL_DECOMPOSITION_H */

