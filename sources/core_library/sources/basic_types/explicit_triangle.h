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

#ifndef EXPLICIT_TRIANGLE_H
#define EXPLICIT_TRIANGLE_H

#include "vertex.h"

/**
 * @brief A class represeting an explicit triangle.
 * The vertices are represented explicitly, i.e., the coordinates of each vertex in the boundary of the triangle are encoded.
 */
class Explicit_Triangle
{
public:
    /**
     * @brief A constructor method
     *
     */
    Explicit_Triangle(){}

    /**
     * @brief A public method that returns the vertex at position pos
     *
     * @param pos an integer representing the position of the vertex in the triangle
     * @return the Vertex at the given position
     */
    inline Vertex& get_vertex(int pos) { return vertices[pos]; }
    /**
     * @brief Public methods that add a vertex to the triangle
     *
     * @param v a Vertex to add to the triangle
     */
    inline void add_vertex(Vertex &v) { vertices.push_back(v); }
    /**
     * @brief Public method that returns the number of vertices of the triangle
     *
     * @return an integer represeting the number of vertices
     */
    inline int vertices_num() const { return vertices.size(); }
    /**
     * @brief
     *
     * @param out
     * @param p
     * @return std::ostream &operator
     */
    inline friend std::ostream& operator<<(std::ostream& out, const Explicit_Triangle& p)
    {
        for(int i=0; i<p.vertices_num(); i++)
        {
            out << p.vertices[i] << endl;
        }
        return out;
    }

private:
    /// A private variable representing the array of vertices forming the triangle
    vector<Vertex> vertices;
};

#endif // EXPLICIT_TRIANGLE_H
