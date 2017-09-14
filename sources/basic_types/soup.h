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

#ifndef SIMPLEXES_SOUP_H
#define SIMPLEXES_SOUP_H

#include "box.h"
#include "triangle.h"
#include "mesh.h"
#include "explicit_triangle.h"

using namespace std;

/**
 * @brief A class representing a soup of triangles in which there is only a list of triangle, representing explcitly the vertices in their boundary
 *
 */
class Soup {
public:
    /**
     * @brief A constructor method
     *
     */
    Soup()
    {
        domain = Box();
    }
    /**
     * @brief A copy-constructor method
     *
     * @param orig
     */
    Soup(const Soup& orig)
    {
        this->domain = orig.domain;
        this->triangles = orig.triangles;
    }
    /**
     * @brief A destructor method
     *
     */
    virtual ~Soup()
    {
        triangles.clear();
    }
    /**
     * @brief A public method that clears the triangles array
     *
     */
    inline void clear() { triangles.clear(); }
    /**
     * @brief A public method that returns the soup domain
     *
     * @return Box&, the mesh domain
     */
    inline Box& get_domain() { return this->domain; }
    /**
     * @brief A public method that sets the soup domain
     *
     * @param d a Box& argument, representing the domain to set
     */
    inline void set_domain(Box& d) { this->domain = d; }
    /**
     * @brief A public method that initializes the triangles array
     *
     * @param num_t a itype argument, representing the number of all triangles
     */
    inline void init_triangles_array(itype num_t) { triangles.reserve(num_t); }
    /**
     * @brief A public method that adds a triangle to the corresponding array
     *
     * @param pos an integer argument representing the dimension of the triangle
     * @param t an ExplicitCell& argument containing the triangle to insert
     */
    inline void add_triangle(Explicit_Triangle& t) { this->triangles.push_back(t); }
    /**
     * @brief A public method that returns a triangle at a given position
     *     
     * @param id an itype argument representing the position index of the triangle in the corresponding array
     * @return ExplicitCell& representing triangle at the given position
     */
    inline Explicit_Triangle& get_triangle(itype id) { return this->triangles[id-1]; }
    /**
     * @brief A public methot that returns the number of triangles encoded by the soup
     *     
     * @return a itype containing the number of triangles
     */
    inline itype get_triangles_num() { return this->triangles.size(); }

protected:
    ///A private varible representing the soup domain
    Box domain;
    ///A private varible representing the top simplexes list of the soup
    vector<Explicit_Triangle> triangles;
};

#endif // SIMPLEXES_SOUP_H
