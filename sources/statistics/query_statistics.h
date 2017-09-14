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

#ifndef _QUERYSTATISTICS_H
#define	_QUERYSTATISTICS_H

#include <list>
#include <vector>
#include <set>
#include <boost/dynamic_bitset.hpp>
#include "basic_types/basic_wrappers.h"
#include "utilities/sorting.h"

using namespace std;
///A class representing a container used to store the statistics obtained from a query over a spatial index
class QueryStatistics{
public:
    ///A public variable representing the number of nodes visited during query
    int numNode;
    ///A public variable representing the number of leaves visited during query
    int numLeaf;
    ///A public variable representing the number of geometric tests executed during a query
    utype numGeometricTest;
    ///A public variable representing the triangles (without duplicates) satisfying a query
    ivect triangles;

    ///A constructor method
    QueryStatistics()  { numNode=numLeaf=numGeometricTest=0; } // used for point locations
    ///A constructor method
    QueryStatistics(itype tri, int perc_res) // used for box queries
    {
        numNode=numLeaf=numGeometricTest=0;
        access_per_tri = ivect(tri,0);

        int reserving = tri / perc_res;
        triangles.reserve(reserving);

        checked_tri = boost::dynamic_bitset<>(tri);
    }
    ///A destructor method
    virtual ~QueryStatistics()
    {
        numNode=numLeaf=numGeometricTest=0;
        checked_tri.reset();
        triangles.clear();
        access_per_tri.clear();
    }

    /**
     * @brief A public procedure that resets the variables of this class (box and line queries wrapper)
     */
    inline void reset(bool )
    {
        numNode=numLeaf=numGeometricTest=0;
        checked_tri.reset();
        triangles.clear();
        fill(access_per_tri.begin(),access_per_tri.end(),0);
    }
    /**
     * @brief A public procedure that resets the variables of this class (point queries wrapper)
     */
    inline void reset()
    {
        numNode=numLeaf=numGeometricTest=0;
        triangles.clear();
    }

    /**
     * @brief is_checked
     * @param t_id
     * @return
     */
    inline bool is_checked(int t_id) { return checked_tri[t_id-1]; }
    /**
     * @brief set_checked
     * @param t_id
     */
    inline void set_checked(int t_id) { checked_tri.set(t_id-1); }
    /**
     * @brief increase_tri_counter
     * @param t_id
     */
    inline void increase_tri_counter(int t_id) { access_per_tri[t_id-1]++; }
    /**
     * @brief get_access_per_tri_array
     * @return
     */
    inline ivect& get_access_per_tri_array() { return access_per_tri; }

private:
    /// I made these variables private in order to provide useful frontends..

    ///A public variable representing the list of triangles that has been checked during a box query
    boost::dynamic_bitset<> checked_tri;
    ///A public array, with an entry for each triangle, containing the number of accesses a triangle had during a query
    ivect access_per_tri;
};

#endif	/* _QUERYSTATISTICS_H */

