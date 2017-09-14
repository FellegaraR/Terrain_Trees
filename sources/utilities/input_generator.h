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

#ifndef _INPUT_GENERATOR_H
#define	_INPUT_GENERATOR_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include "basic_types/point.h"
#include "basic_types/box.h"
#include "geometry/geometry_wrapper.h"

using namespace std;
///A static class that provides an interface for generating random input files for point location, box query and line queries
class Input_Generator {
public:
    ///A public static method that generates a list of points and save it in a file
    /*!
     * \param region a Box& argument, represents the mesh domain in which generate the points
     * \param num_entries an utype representing the number of points to generate
     * \param output a string representing the file name
     */
    static void generate_random_point_inputs(Box &region, utype num_entries, string output);
    ///A public static method that generates a list of points and save it in a set
    /*!
     * \param region a Box& argument, represents the mesh domain in which generate the points
     * \param num_entries an utype representing the number of points to generate
     * \param output a set that will contain the generated points
     */
    static void generate_random_point_inputs(Box &region, utype num_entries, set<Point> &output);
    /**
     * @brief A public static method that generates a list of points near to the centroids of randomly chosen triangle
     *
     * @param region a Box& argument, represents the mesh domain in which generate the points
     * @param num_entries an utype representing the number of points to generate
     * @param mesh a Mesh& argument, representing the triangle mesh
     * @param output a string representing the file name
     */
    static void generate_near_point_inputs(Box& region, utype num_entries, Mesh &mesh, string output);
    ///A public static method that generates a list of boxes and save it in a file
    /*!
     * @param region a Box& argument, represents the mesh domain in which generate inputs
     * @param ratio a coord_type argument, represents the side of the boxes to generate, picked as a percentage of the mesh domain diagonal
     * @param num_entries an utype representing the number of points to generate
     * @param output a string representing the file name
     */
    static void generate_random_box_inputs(Box& region, coord_type ratio, utype num_entries, string output);
    ///A public static method that generates a list of boxes and save it in a set
    /*!
     * @param region a Box& argument, represents the mesh domain in which generate inputs
     * @param ratio a coord_type argument, represents the side of the boxes to generate, picked as a percentage of the mesh domain diagonal
     * @param num_entries an utype representing the number of points to generate
     * \param output a set that will contain the generated boxes
     */
    static void generate_random_box_inputs(Box& region, coord_type ratio, utype num_entries, set<Box> &output);
    /**
     * @brief A public static method that generates a list of boxes with the minimum point of each box corresponding to the centroid of a randomly chosen triangle
     *
     * @param region a Box& argument, represents the mesh domain in which generate inputs
     * @param ratio a coord_type argument, represents the side of the boxes to generate, picked as a percentage of the mesh domain diagonal
     * @param num_entries an utype representing the number of points to generate
     * @param mesh a Mesh& argument, representing the triangle mesh
     * @param output a string representing the file name
     */
    static void generate_near_box_inputs(Box& region, coord_type ratio, utype num_entries, Mesh &mesh, string output);
    ///A public static method that generates a list of lines
    /*!
     * @param region a Box& argument, represents the mesh domain in which generate inputs
     * @param ratio a coord_type argument, represents the lenght of the lines to generate, picked as a percentage of the mesh domain diagonal
     * @param num_entries an utype representing the number of points to generate
     * @param output a string representing the file name
     */
    static void generate_random_line_inputs(Box& region, coord_type ratio, utype num_entries, string output);
    /**
     * @brief A public static method that generates a list of lines with the minimum point of each line corresponding to the centroid of a randomly chosen triangle
     * @param region a Box& argument, represents the mesh domain in which generate inputs
     * @param ratio a coord_type argument, represents the lenght of the lines to generate, picked as a percentage of the mesh domain diagonal
     * @param num_entries an utype representing the number of points to generate
     * @param mesh a Mesh& argument, representing the triangle mesh
     * @param output a string representing the file name
     */
    static void generate_near_line_inputs(Box& region, coord_type ratio, utype num_entries, Mesh &mesh, string output);
private:
    ///A constructor method
    Input_Generator() {}
    ///A copy-constructor method
    Input_Generator(const Input_Generator&) {}
    ///A destructor method
    virtual ~Input_Generator() {}

    static void generate_random_point(Point &p, Box& region);
    static void generate_random_versor(Point &p);

    static void generate_random_boxes(Box &region, set<Box> &boxes, utype num_entries, coord_type edge);
    static void generate_near_boxes(Box &region, set<Box> &boxes, utype num_entries, coord_type edge, Mesh& mesh);

    static void generate_random_lines(Box &region, set<Box> &boxes, utype num_entries, coord_type edge);
    static void generate_near_lines(Box &region, set<Box> &boxes, utype num_entries, coord_type edge, Mesh& mesh);
};

#endif	/* _INPUTGENERATOR_H */

