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

#ifndef _BOX_H
#define	_BOX_H

#include "point.h"
#include "triangle.h"

#include <iostream>

/// A class representing a box in a 3D space
class Box {
public:
    ///A constructor method
    Box()
    {
        min = Point();
        max = Point();
    }
    ///A copy-constructor method
    Box(const Box& orig)
    {
        min = orig.min;
        max = orig.max;
    }
    ///A constructor method
    Box(Point& min, Point& max)
    {
        this->min = min;
        this->max = max;
    }
    ///A destructor method
    virtual ~Box() {}

    ///Public method that returns the minimum point of the box.
    /*!
        \return the reference to the minimum point of the box
    */
    inline Point& get_min() { return min; }
    ///Public method that returns the maximum point of the box.
    /*!
        \return the reference to the maximum point of the box
    */
    inline Point& get_max() { return max; }
    ///Public method that sets the minimum point of the box
    /*!
     * \param x a coord_type argument, representing the x-coordinate
     * \param y a coord_type argument, representing the y-coordinate
    */
    inline void set_min(coord_type x, coord_type y) { this->min.set(x,y); }
    ///Public method that sets the maximum point of the box
    /*!
     * \param x a coord_type argument, representing the x-coordinate
     * \param y a coord_type argument, representing the y-coordinate
    */
    inline void set_max(coord_type x, coord_type y) { this->max.set(x,y); }

    ///Public method that sets the minimum point of the box
    /*!
     * \param p is the point
    */
    inline void set_min(Point& p) { this->min.set(p); }
    ///Public method that sets the maximum point of the box
    /*!
     * \param p is the point
    */
    inline void set_max(Point& p) { this->max.set(p); }
    /**
     * @brief A public method that returns the diagonal of the box
     * @return a coord_type value
     */
    inline coord_type get_diagonal()
    {
        coord_type xedge = fabs(max.get_x()-min.get_x());
        coord_type yedge = fabs(max.get_y()-min.get_y());
        return sqrt(xedge*xedge+yedge*yedge);
    }

    ///Public method that checks if two boxes intersect each other.
    /*!
     * \param other a Box& argument, represents the box we need to test
     * \return a boolean value, true if other intersects the current box, false otherwise
     */
    inline bool intersects(Box& other)
    {
        if (this->max.get_x() < other.min.get_x()) return false;
        if (this->min.get_x() > other.max.get_x()) return false;
        if (this->max.get_y() < other.min.get_y()) return false;
        if (this->min.get_y() > other.max.get_y()) return false;
        return true;
    }
    /**
     * @brief A public procedure that checks if the current box contains completely another box
     * NOTA: the function considers all the faces of the box as closed.
     *
     * @param other a Box& argument, represents the box we need to test
     * @return a boolean value, true if other is completely contained by the current box, false otherwise
     */
    inline bool completely_contains(Box& other)
    {
        if(this->min.get_x() < other.min.get_x() &&
                this->min.get_y() < other.min.get_y() &&                
                this->max.get_x() > other.max.get_x() &&
                this->max.get_y() > other.max.get_y() )
            return true;
        else
            return false;
    }
    ///Public method that checks if a point is inside the current box.
    /*!
     * This method is not used to check if a point is inside the current box, beacause it considers
     * all the six faces of the box as closed.
     *
     * \param p a Point& argument, represents the point to check
     * \return a boolean value, true if the point is inside domain, false otherwise
     */
    inline bool contains_with_all_closed_faces(Point& p)
    {
        for(int i=0; i<p.get_dimension(); i++)
        {
            if(!is_in_range_all_closed(min.get_c(i),max.get_c(i),p.get_c(i)))
                return false;
        }
        return true;
    }
    ///Public method that checks if a point is inside the current box.
    /*!
     * This method considers as closed only three faces, i.e., the ones incidentes in the minimum point.
     * As exceptions, for those boxes (partially) located at the border of the mesh domain,
     * it consider the faces incident in the border as closed.
     *
     * \param p a Point& argument, represents the point to check
     * \param max a Point& argument, represent the maximum point of the mesh domain
     * \return a boolean value, true if the point is inside the box, false otherwise
     */
    inline bool contains(const Point &p, const Point &max)
    {
        //we consider border cases
        // if one of the coordinates are equal to the the one of the mesh domain
        // then we consider that face as closed, otherwise we consider it open
        for(int i=0; i<p.get_dimension(); i++)
        {
            if(!is_in_range(this->min.get_c(i),this->max.get_c(i),p.get_c(i),max.get_c(i)))
                return false;
        }
        return true;
    }
    ///A public method that checks if a point is inside the mesh domain, if the domain doesn't contain the point enlarges it.
    /*!
     * \param p a Point& argument, represents the point to check
     */
    inline void resize(Point& p)
    {
        // If we are big enough to contain the point, do nothing;
        if (contains_with_all_closed_faces(p)) return;

        for(int i=0; i<p.get_dimension(); i++)
        {
            if (p.get_c(i) < this->min.get_c(i)) this->min.set_c(i,p.get_c(i));
            if (p.get_c(i) > this->max.get_c(i)) this->max.set_c(i,p.get_c(i));
        }
    }
    /**
     * @brief operator ==
     * @param p
     * @param q
     * @return
     */
    inline friend bool operator== (const Box& p, const Box& q) { return ((p.min == q.min) && (p.max == q.max)); }
    /**
     * @brief operator !=
     * @param p
     * @param q
     * @return
     */
    inline friend bool operator!= (const Box& p, const Box& q) { return !(p == q); }
    /**
     * @brief operator <
     * @param s
     * @return
     */
    inline bool operator<(const Box& s) const { return ((this->min < s.min) || (this->min == s.min && this->max < s.max)); }
    /**
     * @brief operator >
     * @param s
     * @return
     */
    inline bool operator>(const Box& s) const { return ((this->min > s.min) || (this->min == s.min && this->max > s.max)); }
    /**
     * @brief operator <<
     * @param out
     * @param p
     * @return
     */
    inline friend std::ostream& operator<<(std::ostream& out, const Box& p)
    {
        out << p.min << " " << p.max;
        return out;
    }

    /**
     * @brief A public procedure computing the minimum distance between a point and the box
     * @param p
     * @return
     */
    inline coord_type min_distance(Point &p)
    {
        coord_type sqrDist = 0;

        for(int i=0; i<this->min.get_dimension(); i++)
        {
            coord_type dist;
            if(p.get_c(i) < this->min.get_c(i))
            {
                dist = p.get_c(i) - this->min.get_c(i);
                sqrDist += dist * dist;
            }
            else if(p.get_c(i) > this->max.get_c(i))
            {
                dist = p.get_c(i) - this->max.get_c(i);
                sqrDist += dist * dist;
            }
        }

        return sqrt(sqrDist);
    }

private:
    /// A protected variable representing the minimum point of the box
    Point min;
    /// A protected variable representing the maximum point of the box
    Point max;

    /**
     * @brief A private procedure that checks if a coordinate is in a range
     * NOTA: if the max of the range is equal to the max value of the mesh domain
     * the range is considered as closed, otherwise as open.
     *
     * @param min a coord_type representing the lowest value of the range
     * @param max a coord_type representing the higher value of the range
     * @param coord a coord_type representing the coordinate value
     * @param max_dom a coord_type representing the higher value of the mesh domain
     * @return true if coord is in the range, false otherwise
     */
    inline bool is_in_range(coord_type min, coord_type max, coord_type coord, coord_type max_dom)
    {
        if(max==max_dom){
            if(max < coord){
                return false;
            }
        }
        else if(max <= coord){
            return false;
        }

        if (min > coord){
            return false;
        }

        return true;
    }
    /**
     * @brief A private procedure that checks if a coordinate is in a range
     * NOTA: the range min-max is considered as closed on both sides
     *
     * @param min a coord_type representing the lowest value of the range
     * @param max a coord_type representing the higher value of the range
     * @param coord a coord_type representing the coordinate value
     * @return true if coord is in the range, false otherwise
     */
    inline bool is_in_range_all_closed(coord_type min, coord_type max, coord_type coord)
    {
        if ( (max < coord) || (min > coord) ) return false;
        else return true;
    }


};

#endif	/* _BOX_H */

