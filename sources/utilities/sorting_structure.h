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

#ifndef SORTING_STRUCTURE_H
#define SORTING_STRUCTURE_H

#include <algorithm>
#include <iostream>
#include "basic_types/basic_wrappers.h"

///A container used to store couple of vertex and triangle indexes
class vertex_triangle_pair
{
private:
    ///The vertex index
    itype v;
    ///The triangle index
    itype t;

public:
    vertex_triangle_pair() { v = t = 0; }
    vertex_triangle_pair(itype v, itype t)
    {
        this->v = v;
        this->t = t;
    }

    inline itype get_v() { return this->v; }
    inline itype get_t() { return this->t; }

    bool operator < (const vertex_triangle_pair& p) const { return (v < p.v); }
};

///A container used to store triple of vertex, vertex and triangle indexes
class edge_triangle_tuple
{
private:
    ///The first vertex index
    itype v1;
    ///The second vertex index
    itype v2;
    ///The triangle index
    itype t;
    /// face position in t
    short f_pos;

public:
    edge_triangle_tuple() { v1 = v2 = t = f_pos = 0; }
    edge_triangle_tuple(const edge_triangle_tuple &orig)
    {
        this->v1 = orig.v1;
        this->v2 = orig.v2;
        this->t = orig.t;
        this->f_pos = orig.f_pos;
    }

    inline itype get_v1() { return this->v1; }
    inline itype get_v2() { return this->v2; }
    inline itype get_t() { return this->t; }
    inline short get_f_pos() { return this->f_pos; }
    inline void print() { std::cout<<v1<<" "<<v2<<" "<<t<<" "<<f_pos<<endl; }

    bool operator < (const edge_triangle_tuple& p) const { return (v1 < p.v1) ||
                (v1 == p.v1 && v2 < p.v2) ||
                (v1 == p.v1 && v2 == p.v2 && t < p.t) ; }
    inline bool operator==(const edge_triangle_tuple& s) const { return (v1 == s.v1 && v2 == s.v2); }
    inline bool operator!=(const edge_triangle_tuple& s) const { return !(*this==s); }
    //
    inline void sort_and_set(itype vid1, itype vid2, itype tid)
    {
        v1 = min(vid1,vid2);
        v2 = max(vid1,vid2);
        t = tid;
    }
    inline void sort_and_set(itype vid1, itype vid2, itype tid, short f_p)
    {
        sort_and_set(vid1, vid2, tid);
        f_pos = f_p;
    }
    inline bool has_not(itype v_ind) { return (v_ind != v1 && v_ind != v2); }
};

#endif // SORTING_STRUCTURE_H
