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

#ifndef PRIORITY_QUEUE_ELEMENT_H
#define PRIORITY_QUEUE_ELEMENT_H

#include <queue>
#include "basic_types/box.h"

template<class N> class PQueue_Element
{
public:
    PQueue_Element(coord_type d, N *n, int level, Box &nd) /// constructor method for a tree block
    {
        this->distance = d;
        this->n = n;
        this->value = level;
        this->n_dom = nd;
    }
    PQueue_Element(coord_type d, int id) /// constructor method for a vertex
    {
        this->distance = d;
        this->n = NULL;
        this->value = id;
    }

    inline friend std::ostream& operator<<(std::ostream& out, const PQueue_Element& p)
    {
        out << p.distance << " ";
        if(p.is_block())
            out << *(p.n) << " ";
        else
            out << "v_id: " << p.value << " ";
        return out;
    }

    inline bool is_block() const { return this->n != NULL; }

    inline void set_distance(coord_type d) { this->distance = d; }
    inline void set_domain(Box &d) { this->n_dom = d; }
    inline void set_node(N* n) { this->n = n; }
    inline void set_value(int v) { this->value = v; }

    inline coord_type get_distance() { return this->distance; }
    inline N* get_node() { return n; }
    inline int get_value() { return value; }
    inline Box& get_domain() { return n_dom; }

private:
    N *n;
    Box n_dom;
    int value; // level of the block, or vertex position index
    coord_type distance;
};

struct sort_Element
{
    template<class N> bool operator()(PQueue_Element<N> &s1, PQueue_Element<N> &s2)
    {
        return (s1.get_distance() > s2.get_distance());
    }
};

template<class N> using priority_iNN_queue = priority_queue<PQueue_Element<N>, std::vector<PQueue_Element<N> >, sort_Element>;

#endif // PRIORITY_QUEUE_ELEMENT_H
