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

#include "reader.h"
#include <sstream>
#include <algorithm>
#include "geometry/geometry_wrapper.h"
#include "utilities/string_management.h"

bool Reader::read_mesh(Mesh &mesh, string path)
{
    string extension = string_management::get_file_extension(path);
    if(extension == "tri")
        return Reader::read_mesh_tri(mesh,path);
    else if(extension == "off")
        return Reader::read_mesh_off(mesh,path);
    else
    {
        cerr << "[ERROR] unsopported file format. " << endl;
        return false;
    }
}

bool Reader::read_mesh_off(Mesh& mesh, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    string l;
    getline(input,l); // trow away the first line
    getline(input,l);
    vector<string> lt;
    string_management::tokenize(l,lt," ");
    itype num_vertices = atol(lt.at(0).c_str());
    itype num_triangles = atol(lt.at(1).c_str());

    if (num_vertices == 0 || num_triangles == 0)
    {
        cerr << "This is not a valid .off file: " << path << endl;
        return false;
    }

    mesh.reserve(num_vertices, num_triangles);

    Reader::read_vertices_list(mesh,input,num_vertices);

    itype v[3], index;
    int num_v;
    for (int i = 0; i < num_triangles; i++)
    {
        input >> num_v;

        if(num_v != 3)
        {
            cerr << "[ERROR] the input mesh must be a pure triangle mesh. read a simplex with "<< num_v << "vertices." << endl;
            return false;
        }

        for (int j = 0; j < num_v; j++)
        {
            input >> index;
            v[j] = index+1;
        }
        Triangle t = Triangle(v[0], v[1], v[2]);
        mesh.add_triangle(t);
    }

    return true;
}

bool Reader::read_mesh_tri(Mesh& mesh, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    string l;
    getline(input,l);
    vector<string> lt;
    string_management::tokenize(l,lt," ");
    itype num_vertices = atoi(lt.at(0).c_str());

    if (num_vertices == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    mesh.reserve_vertices_space(num_vertices);

    Reader::read_vertices_list(mesh,input,num_vertices);

    itype num_triangles;
    input >> num_triangles;

    if (num_triangles == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    mesh.reserve_triangles_space(num_triangles);

    itype v[3], index;
    for (itype i = 0; i < num_triangles; i++)
    {
        for (int j = 0; j < 3; j++) {
            input >> index;
            v[j] = index+1;
        }
        Triangle t = Triangle(v[0], v[1], v[2]);
        mesh.add_triangle(t);
    }

    return true;
}

bool Reader::read_vertices(Mesh& mesh, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    string extension = string_management::get_file_extension(path);
    itype num_vertices;
    vector<string> lt;
    string l;

    getline(input,l);

    if(extension == "tri")
    {
        string_management::tokenize(l,lt," ");
        num_vertices = atoi(lt.at(0).c_str());

        if (num_vertices == 0)
        {
            cerr << "This is not a valid .tri file: " << path << endl;
            return false;
        }
    }
    else if(extension == "off")
    {
         getline(input,l);
         string_management::tokenize(l,lt," ");
         num_vertices = atoi(lt.at(0).c_str());

        if (num_vertices == 0)
        {
            cerr << "This is not a valid .off file: " << path << endl;
            return false;
        }
    }

    mesh.reserve_vertices_space(num_vertices);
    Reader::read_vertices_list(mesh,input,num_vertices);

    return true;
}

void Reader::read_vertices_list(Mesh &mesh, ifstream &input, itype num_vertices)
{
    bool is2D = false, isMulti = false;
    string line;
    //read the vertices and set the mesh domain
    for (itype i = 0; i < num_vertices; i++)
    {
        getline(input,line);
        vector<string> line_tokens;
        string_management::tokenize(line,line_tokens," ");

        if (input.eof())
            break;

        Vertex v;
        if(line_tokens.size() == 3) // 3D
            v = Vertex(atof(line_tokens[0].c_str()),atof(line_tokens[1].c_str()),atof(line_tokens[2].c_str()));
        else if(line_tokens.size() == 2) // 2D
        {
            v = Vertex(atof(line_tokens[0].c_str()),atof(line_tokens[1].c_str()));
            is2D = true;
        }
        else if(line_tokens.size() > 3) // MULTIFIELD
        {
            isMulti=true;
            dvect fields;
            for(int f=2; f<line_tokens.size(); f++)
                fields.push_back(atof(line_tokens[f].c_str()));
            v = Vertex(atof(line_tokens[0].c_str()),atof(line_tokens[1].c_str()),fields);
        }
        mesh.add_vertex(v);
        if (i == 0) {
            Box b = Box(v, v);
            mesh.set_domain(b);
        } else {
            mesh.get_domain().resize(v);
        }
    }

    if(is2D)
        cerr<<"[NOTA] The points are embedded in a 2D space."<<endl;
    else// if(!isMulti)
        cerr<<"[NOTA] The points are embedded in a 3D space."<<endl;
    if(isMulti)
        cerr<<"[NOTA] Multifield triangle mesh."<<endl;
}

bool Reader::read_soup(Soup& soup, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    bool is_first = true;

    string line;
    while(input.good())
    {
        Explicit_Triangle tri;
        getline(input,line);        
        if (input.eof())
            break;
        vector<string> line_tokens;
        string_management::tokenize(line,line_tokens,"\t");
        for(unsigned i=0; i<line_tokens.size(); i++)
        {
            vector<string> point_tokens;
            string_management::tokenize(line_tokens[i],point_tokens,",");
            Vertex v;
            if(point_tokens.size() == 3)
                 v = Vertex(atof(point_tokens[0].c_str()),atof(point_tokens[1].c_str()),atof(point_tokens[2].c_str()));
            else if(point_tokens.size() == 2)
                v = Vertex(atof(point_tokens[0].c_str()),atof(point_tokens[1].c_str()));
            tri.add_vertex(v);
            if (is_first)
            {
                Box b = Box(v, v);
                soup.set_domain(b);
                is_first = false;
            } else
            {
                soup.get_domain().resize(v);
            }
        }
        soup.add_triangle(tri);
    }

    input.close();

    return true;
}

void Reader::read_queries(vector<Point>& points, string fileName)
{
    ifstream input(fileName.c_str());
    int size = 0;
    input >> size;
    points.reserve(size);
    coord_type x = 0, y = 0;
    while (input)
    {
        input >> x;
        input >> y;
        if (input.eof())
            break;
        points.push_back(Point(x, y));
    }
}

void Reader::read_queries(vector<Box> &boxes, string fileName)
{
    ifstream input(fileName.c_str());
    int size = 0;
    input >> size;
    boxes.reserve(size);
    coord_type x1 = 0, y1 = 0, x2 = 0, y2 = 0;
    while (input) {
        input >> x1;
        input >> y1;
        input >> x2;
        input >> y2;
        if (input.eof())
            break;
        Point min = Point(x1, y1);
        Point max = Point(x2, y2);
        Box b = Box(min, max);
        boxes.push_back(b);
    }
}

void Reader::read_leaf(Node_T* n, ifstream& input, vector<string>& tokens)
{
    string line;
    int numTri = atoi(tokens.at(1).c_str());
    if(numTri > 0)
    {
        getline(input, line);
        vector<string> tokens2;
        istringstream iss2(line);
        copy(istream_iterator<string > (iss2),
             istream_iterator<string > (),
             back_inserter<vector<string> >(tokens2));
        for (unsigned int i = 1; i < tokens2.size(); i++)
            n->add_triangle(atol(tokens2.at(i).c_str()));
    }
}

void Reader::read_leaf(Node_V *n, ifstream &input, vector<string>& tokens)
{
    string line;
    int numVertex = atoi(tokens.at(1).c_str());
    int numTetra = atoi(tokens.at(2).c_str());

    if (numVertex > 0)
    {
        getline(input, line);
        vector<string> tokens2;
        istringstream iss2(line);
        copy(istream_iterator<string > (iss2),
             istream_iterator<string > (),
             back_inserter<vector<string> >(tokens2));
        for (unsigned int i = 1; i < tokens2.size(); i++)
            n->add_vertex(atol(tokens2.at(i).c_str()));
    }

    if(numTetra > 0)
    {
        getline(input, line);
        vector<string> tokens3;
        istringstream iss3(line);
        copy(istream_iterator<string > (iss3),
             istream_iterator<string > (),
             back_inserter<vector<string> >(tokens3));
        for (unsigned int i = 1; i < tokens3.size(); i++)
            n->add_triangle(atol(tokens3.at(i).c_str()));
    }
}

bool Reader::read_noised_field_value(Mesh& mesh, string mesh_name, dvect &original_vertex_fields)
{
    original_vertex_fields.assign(mesh.get_vertices_num(),0);
    stringstream ss;
    ss << mesh_name << ".field";
    ifstream input(ss.str().c_str());
    coord_type f;

    if (input.is_open() == false) {
        cerr << "The mesh is already noised." << endl;
        return false;
    }

    for(int v=1; v<=mesh.get_vertices_num(); v++)
    {
        Vertex& vert = mesh.get_vertex(v);
        original_vertex_fields[v-1] = vert.get_z();
        input >> f;
        vert.set_c(2,f);
    }

    return true;
}
