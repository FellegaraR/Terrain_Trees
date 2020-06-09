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

#ifndef _VERTEX_H
#define	_VERTEX_H

#include "point.h"

using namespace std;

///An inner-class, extending Point, representing a vertex in a triangle mesh
class Vertex : public Point
{
public:
    ///A constructor method
    Vertex() : Point() { this->fields.push_back(0); }
    ///A copy-constructor method
    Vertex(const Vertex& orig) : Point(orig) { this->fields=orig.fields; }
    ///A constructor method
    /*!
     * \param x a coord_type argument, representing the x coordinate
     * \param y a coord_type argument, representing the y coordinate
     * \param field a coord_type argument, representing the vertex field
     */
    Vertex(coord_type x, coord_type y, coord_type field) : Point(x,y)
    {
       this->fields.push_back(field);
    }
    ///A constructor method
    /*!
     * \param x a coord_type argument, representing the x coordinate
     * \param y a coord_type argument, representing the y coordinate
     */
    Vertex(coord_type x, coord_type y) : Point(x,y)
    {
       this->fields.push_back(0);
    }
    ///A constructor method
    /*!
     * \param x a coord_type argument, representing the x coordinate
     * \param y a coord_type argument, representing the y coordinate
     * \param fields a vector argument, representing the vertex fields
     */
    Vertex(coord_type x, coord_type y, dvect &fields) : Point(x,y)
    {
       this->fields = fields;
    }
    ///A destructor method
    virtual ~Vertex() { }
    ///A public method that returns the vertex elevation
    /*!
     * \return an integer, representing the elevation
     */
    inline coord_type get_z() { return fields[0]; }
    /**
     * @brief A public method that returns one of the coordinates or the elevation value of the vertex
     *
     * @param pos an integer representing the coordinate position in the point array
     * @return a coord_type representing the coordinate
     */
    inline coord_type get_c(int pos) const
    {
        if(pos == this->get_dimension())
            return this->fields[0];
        else
            return this->coords[pos];
    }
    /**
     * @brief A public method that set the value of a coordinate or the elevation value
     *
     * @param pos an integer representing the coordinate position in the point array
     * @param c a coord_type representing the coordinate value
     */
    inline void set_c(int pos, coord_type c)
    {
        if(pos == this->get_dimension())
            this->fields[0] = c;
        else
            this->coords[pos] = c;
    }
    /**
     * @brief operator <<
     * @param out
     * @param p
     * @return
     */    
    inline friend std::ostream& operator<<(std::ostream& out, const Vertex& p)
    {
        out << p.coords[0] << " " << p.coords[1]<<" ";
        for(auto f : p.fields)
            out << f << " ";
        return out;
    }

    /**
     * @brief A public procedure that computes the norm value of the vector this-v
     *
     * @param v a Vertex& representing the other side of the vector
     * @return the norm value
     */
    inline coord_type norm(Vertex& v)
    {
        coord_type xdist = v.get_x()-this->coords[0];
        coord_type ydist = v.get_y()-this->coords[1];
        coord_type zdist = v.get_z()-this->fields[0];

        return sqrt(xdist*xdist + ydist*ydist + zdist*zdist);
    }

    /**
     * @brief A public procedure returning the scalar product of the vectors v1-this and v2-this
     *
     * @param v1 a Vertex& representing a vertex
     * @param v2 a Vertex& representing a vertex
     * @return the scalar product value
     */
    inline coord_type scalar_product(Vertex& v1,Vertex& v2)
    {
        return(((v1.get_x()-this->coords[0])*(v2.get_x()-this->coords[0]))
               +((v1.get_y()-this->coords[1])*(v2.get_y()-this->coords[1]))
               +((v1.get_z()-this->fields[0])*(v2.get_z()-this->fields[0])));
    }

    /**
     * @brief A public procedure that adds a field value in the last position of the fields array
     * @param f a coord_type variable representing the field value to add
     */
    inline void add_field(coord_type f) { fields.push_back(f); }
    /**
     * @brief A public procedure that return the field value stored at a certain position
     * @param pos an integer value, encoding the position in the fields array
     * @return the field value
     */
    inline coord_type get_field(int pos) { return fields[pos]; }
    /**
     * @brief A public procedure returning the number of field values encoded in the vertex
     * @return the number of field value
     */
    inline int get_fields_num() { return fields.size(); }
    
    void set_gradientMatrix(FG gradient){this->gradientMatrix.push_back(gradient);}
    void print_gradientMatrix(){
        for(int i=0;i<this->gradientMatrix.size();i++)
            cout<<"stored gradient of field "<<i<<":"<<gradientMatrix[i][0]<<", "<<gradientMatrix[i][1]<<endl;
    
    };
    
    vect_FG& get_gradientMatrix(){ return this->gradientMatrix;
    
    
    };
private:
    ///A private variable that represents the field values
    /// NOTA: the first field value is always the z-value (or elevation)
    dvect fields;
    vect_FG gradientMatrix;
};

#endif	/* _VERTEX_H */

