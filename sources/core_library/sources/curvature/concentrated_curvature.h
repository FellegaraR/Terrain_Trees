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

#ifndef CONCENTRATED_CURVATURE_H
#define CONCENTRATED_CURVATURE_H

#include "abstract_curvature.h"
#include "border_checker.h"

enum ConcentratedType {CONCENTRATED = 0, GAUSS_ANGLE_DEFICIT = 1};

class Concentrated_Curvature : public Abstract_Curvature, public Border_Checker
{
public:
    //if divide_by_area==1 then gaussian angle deficit
    //if divide_by_area==0 then concentrated curvature
    Concentrated_Curvature(ConcentratedType type) : Abstract_Curvature(), Border_Checker()
    {
        this->type = type;
    }

private:
    //flag if we want to divide by area of the neighborhood or not
    ConcentratedType type;
    leaf_VT vts;
    /**
     * @brief curvature_leaf
     * @param n
     * @param mesh
     */
    void curvature_leaf(Node_V& n, Mesh& mesh);
    /**
     * @brief curvature_leaf
     * @param n
     * @param mesh
     */
    void curvature_leaf(Node_T& n, Box &n_dom, Mesh& mesh);

    void update_local_curvature(itype v_id, itype local_v_id, Triangle &t, itype v_pos,
                                boost::dynamic_bitset<> &is_v_border, dvect &local_curvatures, Mesh &mesh);
    void finalize_local_curvatures(itype v_start, leaf_VT &vts, boost::dynamic_bitset<> &is_v_border, dvect &local_curvatures, Mesh &mesh);
};

#endif // CONCENTRATED_CURVATURE_H
