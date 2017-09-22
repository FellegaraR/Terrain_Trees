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

#ifndef _POINT_H
#define	_POINT_H

#include <iostream>
#include <cmath>

#include "basic_wrappers.h"

using namespace std;

///A class representing a 2D point
class Point {
public:
    ///A constructor method
    Point() { this->coords = { 0, 0 }; }
    ///A copy-constructor method
    Point(const Point& orig) { this->coords = orig.coords; }
    ///A constructor method
    /*!
     * \param x a coord_type argument, representing the x coordinate
     * \param y a coord_type argument, representing the y coordinate
     */
    Point(coord_type x, coord_type y) { this->coords = { x, y }; }
    ///A destructor method
    virtual ~Point() {}
    /**
     * @brief operator ==
     * @param p
     * @param q
     * @return
     */
    inline friend bool operator== (const Point& p, const Point &q) { return ((p.coords[0] == q.coords[0]) && (p.coords[1] == q.coords[1])); }
    /**
     * @brief operator !=
     * @param p
     * @param q
     * @return
     */
    inline friend bool operator!= (const Point& p, const Point &q) { return !(p == q); }
    /**
     * @brief operator <
     * @param s
     * @return
     */
    inline bool operator<(const Point& s) const
    {
        return ((this->coords[0] < s.coords[0]) ||
                (this->coords[0] == s.coords[0] && this->coords[1] < s.coords[1]));
    }
    /**
     * @brief operator >
     * @param s
     * @return
     */
    inline bool operator>(const Point& s) const
    {
        return ((this->coords[0] > s.coords[0]) ||
                (this->coords[0] == s.coords[0] && this->coords[1] > s.coords[1]));
    }
    /**
     * @brief operator +
     * @param s
     * @return
     */
    inline Point operator+(const Point &s) const { return Point(this->coords[0]+s.get_x(),this->coords[1]+s.get_y()); }
    /**
     * @brief operator -
     * @param s
     * @return
     */
    inline Point operator-(const Point &s) const { return Point(this->coords[0]-s.get_x(),this->coords[1]-s.get_y()); }
    /**
     * @brief operator *
     * @param f
     * @return
     */
    inline Point operator*(const coord_type &f) const { return Point(this->coords[0]*f,this->coords[1]*f); }
    ///
    inline friend std::ostream& operator<<(std::ostream& out, const Point& p)
    {
        out << p.coords[0] << " " << p.coords[1];
        return out;
    }
    ///A public method that returns the x coordinate
    /*!
     * \return a coord_type representing the x coordinate
     */
    inline coord_type get_x() const { return this->coords[0]; }
    ///A public method that returns the y coordinate
    /*!
     * \return a coord_type representing the y coordinate
     */
    inline coord_type get_y() const { return this->coords[1]; }
    /**
     * @brief A public procedure returning the coordinate at a given position
     *
     * @param pos an integer representing the coordinate position in the point array
     * @return a coord_type representing the coordinate
     */
    inline coord_type get_c(int pos) const { return this->coords[pos]; }
    /**
     * @brief A public procedure that initializes a coordinate of the point
     *
     * @param pos an integer representing the coordinate position in the point array
     * @param c a coord_type representing the coordinate value
     */
    inline void set_c(int pos, coord_type c) { this->coords[pos] = c; }
    /**
     * @brief A public method that initializes all the coordinates of a point
     * @param x a coord_type representing the value on the x-axis
     * @param y a coord_type representing the value on the y-axis
     */
    inline void set(coord_type x, coord_type y)
    {
        this->coords[0]=x;
        this->coords[1]=y;
    }
    /**
     * @brief A public method that initializes all the coordinates of a point
     * @param p a Point& argument containing the the coordinates value to set
     */
    inline void set(Point &p) { set(p.get_x(),p.get_y()); }
    /**
     * @brief A public method returning the dimension of the point
     * @return an integer representing the dimension
     */
    inline int get_dimension() const { return coords.size(); }

    /**
     * @brief A public method that computes the distance between two points
     *
     * @param v, a Point&
     * @return coord_type, the distance value
     */
    inline coord_type distance(Point& v)
    {
        coord_type distance_factor = 0;
        for(int i=0; i<v.get_dimension(); i++)
        {
            distance_factor += (coords[i]-v.get_c(i))*(coords[i]-v.get_c(i));
        }
        return sqrt(distance_factor);
    }

protected:
    /// A protected array representing the point x and y coordinates
    dvect coords;
};

#endif	/* _POINT_H */

