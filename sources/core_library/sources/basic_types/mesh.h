/*
    This file is part of the triangle Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The trianglesl Trees library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The trianglesl Trees library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the trianglesl Trees library.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _MESH_H
#define	_MESH_H

#include <vector>
#include <climits>

#include "vertex.h"
#include "triangle.h"
#include "box.h"
#include "Matrix.h"

using namespace std;
/// A class representing a triangulated terrain.
/// NOTA: the range of vertices and top cells indexes is defined from 1 to #num_entities (included)
/// 0 is a special value not used to index any element.
class Mesh
{
public:
    ///A constructor method
    Mesh()
    {
        domain = Box();
        vertices = vector<Vertex>();
        triangles = vector<Triangle>();
    }
    ///A copy-constructor method
    Mesh(const Mesh& orig)
    {
        this->triangles = orig.triangles;
        this->vertices = orig.vertices;
        this->domain = orig.domain;
    }
    ///A destructor method
    virtual ~Mesh()
    {
        vertices.clear();
        triangles.clear();
    }
    ///A public method that returns the vertex at the i-th position in the mesh list
    /*!
     * \param id an itype argument, representing the position in the list
     * \return a Vertex&, the vertex at the id-th position in the list
     */
    inline Vertex& get_vertex(itype id) { return this->vertices[id-1]; }

    inline vector<Vertex> get_vertices_array(){return this->vertices;}
    ///A public method that returns the triangle at the i-th position in the mesh list
    /*!
     * \param id an itype argument, representing the position in the list
     * \return a Triangle&, the triangle at the id-th position in the list
     */
    inline Triangle& get_triangle(itype id) { return this->triangles[id-1]; }
    ///A public method that returns the mesh domain
    /*!
     * \return a Box&, the mesh domain
     */
    inline Box& get_domain() { return this->domain; }
    ///A public method that returns the number of mesh vertices
    /*!
     * \return an itype, representing the number of vertices
     */
    inline utype get_vertices_num() { return this->vertices.size(); }
    ///A public method that returns the number of mesh triangles
    /*!
     * \return an itype, representing the number of triangles
     */
    inline utype get_triangles_num() { return this->triangles.size(); }
    ///A public method that sets the mesh domain
    /*!
     * \param d a Box& argument, representing the domain to set
     */
    inline void set_domain(Box& d) { this->domain = d; }
    ///A public method that adds a vertex to the vertices list
    /*!
     * \param v a Vertex& argument, representing the vertex to add
     */
    inline void add_vertex(Vertex& v) { this->vertices.push_back(v); }
    ///A public method that adds a triangle to the triangles list
    /*!
     * \param t a Triangle& argument, representing the triangle to add
     */
    inline void add_triangle(Triangle& t) { this->triangles.push_back(t); }
    ///A public method that initializes the space needed by the vertices and triangles arrays
    /*!
     * \param numV an itype, represents the number of mesh vertices
     * \param numT an itype, represents the number of mesh triangles
     */
    inline void reserve(itype numV, itype numT)
    {
        this->vertices.reserve(numV);
        this->triangles.reserve(numT);
    }
    ///A public method that initializes the space needed by the vertices array
    /*!
     * \param numV an itype, represents the number of mesh vertices
     */
    inline void reserve_vertices_space(itype numV) { this->vertices.reserve(numV); }
    ///A public method that resets the vertices array
    inline void reset_vertices() { this->vertices.clear(); }
    ///A public method that initializes the space needed by the triangles array
    /*!
     * \param numT an itype, represents the number of mesh triangles
     */
    inline void reserve_triangles_space(itype numT) { this->triangles.reserve(numT); }
    ///A public method that resets the triangles array
    inline void reset_triangles() { this->triangles.clear(); }

    /**
     * @brief A public method that return the vertex with the maximum elevation into the triangle
     * @param t a Triangle variable
     * @return the position index of the vertex with the maximum elevation
     */
    inline itype get_max_elevation_vertex(Triangle &t)
    {
        itype max = t.TV(0);
        coord_type max_field = this->get_vertex(t.TV(0)).get_z();
        for(itype v=1; v<t.vertices_num(); v++)
        {
            if(this->get_vertex(t.TV(v)).get_z() > max_field)
            {
                max = t.TV(v);
                max_field = this->get_vertex(t.TV(v)).get_z();
            }
        }
        return max;
    }
    /**
     * @brief A public method that return the vertex with the maximum elevation into the simplex
     * @param vect a vector variable representing the simplex
     * @return the position index of the vertex with the maximum elevation
     */
    inline itype get_max_elevation_vertex(const ivect &vect)
    {
        if(vect.size() == 0)
        {
            cerr << "[ERROR] get_max_field_vertex -> empty vector" << endl;
            return -1;
        }

        itype max = vect[0];
        coord_type max_field = this->get_vertex(vect[0]).get_z();
        for(unsigned v=1; v<vect.size(); v++)
        {
            if(this->get_vertex(vect[v]).get_z() > max_field)
            {
                max = vect[v];
                max_field = this->get_vertex(vect[v]).get_z();
            }
        }
        return max;
    }

    /**
     * @brief A public method that checks if the vertex at position v is flagged as deleted
     *
     * @param v the position of the vertex in the vertices array
     * @return bool, true if the vertex is deleted, false otherwise
     */
    inline bool is_vertex_removed(itype v) { return (get_vertex(v).get_c(0)==INFINITY); }
    /**
     * @brief A public method that checks if the vertex v is flagged as deleted
     *
     * @param v the vertex to check
     * @return bool, true if the vertex is deleted, false otherwise
     */
    inline bool is_vertex_removed(Vertex &v) { return (v.get_c(0)==INFINITY); }
    /**
     * @brief A public method that checks if the top cell of dimension dim at position t is flagged as deleted
     *
     * @param t the top cell position index in that array
     * @return bool, true if the top cell is deleted, false otherwise
     */
    inline bool is_triangle_removed(itype t) { return (get_triangle(t).TV(0) == 0); }
    /**
     * @brief A public method that checks if the top cell t is flagged as deleted
     *
     * @param t the top cell to check
     * @return bool, true if the top cell is deleted, false otherwise
     */
    inline bool is_triangle_removed(Triangle &t) { return (t.TV(0) == 0); }
    /**
     * @brief A public method that flags as removed the vertex with position index v_id
     *
     * @param v_id the position of the vertex to delete
     */
    inline void remove_vertex(itype v_id)
    {
        Vertex& v = get_vertex(v_id);
        v.set_c(0,INFINITY);
    }
    /**
     * @brief A public method that flags as removed the top cell in the dim-th array and with position index t
     *
     * @param dim, the array position
     * @param t the position index of the top cell to delete
     */
    inline void remove_triangle(itype t)
    {
        Triangle& tri = get_triangle(t);
        tri.setTV(0,0);
    }

    /**
     * @brief A public method that swaps the position of two vertices
     * @param v_id1, the position index of the first vertex
     * @param v_id2, the position index of the secton vertex
     */
    inline void vertices_swap(int v_id1, int v_id2) { swap(vertices[v_id1-1],vertices[v_id2-1]); }
    /**
     * @brief A public method that swaps the position of two triangles
     * @param t_id1, the position index of the first triangle
     * @param t_id2, the position index of the second triangle
     */
    inline void triangles_swap(int t_id1, int t_id2) { swap(this->triangles[t_id1-1],this->triangles[t_id2-1]); }
    /**
     * @brief A public method that resizes the top d-cells array.
     * If new_size is smaller than the current size of the array, the position encoded after size are effectively deleted from the array.
     * @param new_size an integer representing the new size of the array
     */
    inline void resize_triangle_array(int new_size) {this->triangles.resize(new_size);}

    inline std::vector<Triangle>::iterator get_t_array_begin(){return this->triangles.begin();}
    inline std::vector<Triangle>::iterator get_t_array_end(){return this->triangles.end();}

    inline void computeInitialQEM(vector<Matrix>* vQEM, vector<dvect >* planes){

  for (int i = 0; i < get_triangles_num(); i++)
    {
        /* faces are triangles */
        for (int j = 0; j < 3; j++)
        {
            double* a = &((*planes)[i][0]);
            (*vQEM)[ get_triangle(i).TV(j) ] += Matrix(a);
        }
    }

    }
    inline void computeTrianglesPlane(vector<dvect >* trPl)
    {
  double coords[3][3];

    for(int i=0; i<get_triangles_num(); i++){

        for(int v=0; v<3; v++){
            coords[0][v] = get_vertex(get_triangle(i).TV(v)).get_x();
            coords[1][v] = get_vertex(get_triangle(i).TV(v)).get_y();
            coords[2][v] = get_vertex(get_triangle(i).TV(v)).get_z();
        }

        double a,b,c,m;

        a = (coords[1][1] - coords[1][0]) * (coords[2][2] - coords[2][0]) - (coords[2][1] - coords[2][0]) * (coords[1][2] - coords[1][0]);

        b = (coords[2][1] - coords[2][0]) * (coords[0][2] - coords[0][0]) - (coords[0][1] - coords[0][0]) * (coords[2][2] - coords[2][0]);

        c = (coords[0][1] - coords[0][0]) * (coords[1][2] - coords[1][0]) - (coords[1][1] - coords[1][0]) * (coords[0][2] - coords[0][0]);

        m = sqrt(a*a + b*b + c*c);
        a = a/m;
        b = b/m;
        c = c/m;

        (*trPl)[i][0]=a;
        (*trPl)[i][1]=b;
        (*trPl)[i][2]=c;
        (*trPl)[i][3]= -1*(a*coords[0][0] + b*coords[1][0] + c*coords[2][0]);
    }

    }

    inline double vertex_error(Matrix q, double x, double y, double z)
{
    return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x + q[5]*y*y
        + 2*q[6]*y*z + 2*q[7]*y + q[10]*z*z + 2*q[11]*z + q[15];
}


    inline double compute_error(int v1, int v2, vector<Matrix>* vQEM, dvect* new_vertex){
      double min_error;
    Matrix q_bar;
    Matrix q_delta;
    assert(new_vertex != NULL);

    /* computer quadric of virtual vertex vf */
    q_bar = (*vQEM)[v1] + (*vQEM)[v2];


    q_delta = Matrix( q_bar[0], q_bar[1],  q_bar[2],  q_bar[3],
                      q_bar[4], q_bar[5],  q_bar[6],  q_bar[7],
                      q_bar[8], q_bar[9], q_bar[10], q_bar[11],
                             0,        0,	      0,        1);

    double vx1 = vertices[v1].get_x();
    double vy1 = vertices[v1].get_y();
    double vz1 = vertices[v1].get_z();

    double vx2 = vertices[v2].get_x();
    double vy2 = vertices[v2].get_y();
    double vz2 = vertices[v2].get_z();


        double vx3 = double (vx1+vx2)/2.0;
        double vy3 = double (vy1+vy2)/2.0;
        double vz3 = double (vz1+vz2)/2.0;

        double error1 = vertex_error(q_bar, vx1, vy1, vz1);
        double error2 = vertex_error(q_bar, vx2, vy2, vz2);
        double error3 = vertex_error(q_bar, vx3, vy3, vz3);

        min_error = std::min(error1, std::min(error2, error3));
        if (error1 == min_error) { (*new_vertex)[0] = vx1; (*new_vertex)[1] = vy1, (*new_vertex)[2] = vz1; }
        else if (error2 == min_error) { (*new_vertex)[0] = vx2; (*new_vertex)[1] = vy2, (*new_vertex)[2] = vz2; }
        else{ (*new_vertex)[0] = vx3; (*new_vertex)[1] = vy3, (*new_vertex)[2] = vz3; }

//    }


    min_error = vertex_error(q_bar, (*new_vertex)[0], (*new_vertex)[1], (*new_vertex)[2]);

    return min_error;
    }
private:
    ///A private varible representing the mesh domain
    Box domain;
    ///A private varible representing the vertices array of the mesh
    vector<Vertex> vertices;
    ///A private varible representing the triangles array of the mesh
    vector<Triangle> triangles;
};

#endif	/* _MESH_H */

