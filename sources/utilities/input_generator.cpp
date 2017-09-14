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

#include "io/writer.h"
#include <cstdlib>
#include <sstream>
#include "input_generator.h"

void Input_Generator::generate_random_point_inputs(Box &region, utype num_entries, string output)
{
    set<Point> points;
    Point p;
    while(points.size()<num_entries)
    {
        Input_Generator::generate_random_point(p,region);
        pair<set<Point>::iterator,bool> ret = points.insert(p);
        if(ret.second)
        {
            cout<<points.size()<<" "<<flush;
        }

    }
    cout<<endl;
    Writer::write_point_queries(points,output + "_point.pqin");
}

void Input_Generator::generate_random_point_inputs(Box &region, utype num_entries, set<Point> &output)
{
    Point p;
    while(output.size()<num_entries)
    {
        Input_Generator::generate_random_point(p,region);
        pair<set<Point>::iterator,bool> ret = output.insert(p);
        if(ret.second)
        {
            cout<<output.size()<<" "<<flush;
        }

    }
    cout<<endl;
}

void Input_Generator::generate_near_point_inputs(Box& region, utype num_entries, Mesh &mesh, string output)
{
    srand((utype)time(NULL));
    set<Point> points;
    iset t_ids;
    Point centroid = Point();
    itype rand_t_id;
    while(points.size()<num_entries)
    {
        rand_t_id = rand() % mesh.get_triangles_num();
        pair<iset::iterator,bool> ret1 = t_ids.insert(rand_t_id);
        if(ret1.second)
        {
            Geometry_Wrapper::get_triangle_centroid(rand_t_id,centroid,mesh);
            if(!region.contains_with_all_closed_faces(centroid))
                continue;
            pair<set<Point>::iterator,bool> ret = points.insert(centroid);
            if(ret.second)
            {
                cout<<points.size()<<" "<<flush;
            }
        }
    }
    cout<<endl;
    Writer::write_point_queries(points,output + "_point.pqin");
}

void Input_Generator::generate_random_point(Point &p, Box& region)
{
    srand(time(NULL));
    coord_type c;
    for(int i=0; i<p.get_dimension(); i++)
    {
        c = (coord_type)(((coord_type)rand())/((coord_type)RAND_MAX)*(region.get_max().get_c(i)-region.get_min().get_c(i))+region.get_min().get_c(i));
        p.set_c(i,c);
    }
}

void Input_Generator::generate_random_versor(Point &p)
{
    srand(time(NULL));
    coord_type c;
    for(int i=0; i<p.get_dimension(); i++)
    {
        c = (coord_type)(((coord_type)rand())/((coord_type)RAND_MAX));
        p.set_c(i,c);
    }
}

void Input_Generator::generate_random_box_inputs(Box& region, coord_type ratio, utype num_entries, string output)
{
    coord_type edge = region.get_diagonal()*ratio;

    set<Box> boxes;
    stringstream ost;
    ost << ratio;

    Input_Generator::generate_random_boxes(region,boxes,num_entries,edge);

    Writer::write_box_queries(boxes,output + "_box_"+ost.str()+".bqin");
}

void Input_Generator::generate_random_box_inputs(Box& region, coord_type ratio, utype num_entries, set<Box> &output)
{
    coord_type edge = region.get_diagonal()*ratio;
    Input_Generator::generate_random_boxes(region,output,num_entries,edge);
}

void Input_Generator::generate_near_box_inputs(Box& region, coord_type ratio, utype num_entries, Mesh& mesh, string output)
{
    coord_type edge = region.get_diagonal()*ratio;

    set<Box> boxes;
    stringstream ost;
    ost << ratio;

    Input_Generator::generate_near_boxes(region,boxes,num_entries,edge,mesh);

    Writer::write_box_queries(boxes,output + "_box_"+ost.str()+".bqin");
}

void Input_Generator::generate_random_line_inputs(Box& region, coord_type ratio, utype num_entries, string output)
{
    coord_type edge = region.get_diagonal()*ratio;

    set<Box> boxes;
    stringstream ost;
    ost << ratio;

    Input_Generator::generate_random_boxes(region,boxes,num_entries,edge);

    Writer::write_box_queries(boxes,output + "_line_"+ost.str()+".lqin");
}

void Input_Generator::generate_near_line_inputs(Box& region, coord_type ratio, utype num_entries, Mesh& mesh, string output)
{
    coord_type edge = region.get_diagonal()*ratio;

    set<Box> boxes;
    stringstream ost;
    ost << ratio;

    Input_Generator::generate_near_boxes(region,boxes,num_entries,edge,mesh);

    Writer::write_box_queries(boxes,output + "_line_"+ost.str()+".lqin");
}

void Input_Generator::generate_random_boxes(Box &region, set<Box> &boxes, utype num_entries, coord_type edge)
{
    Point min;
    Point max;

    while(boxes.size()<num_entries)
    {
        Input_Generator::generate_random_point(min,region);
        max.set((min.get_x() + edge),(min.get_y() + edge));

        //check that max point is into the domain
        if(!region.contains_with_all_closed_faces(max))
            continue;

        Box b = Box(min,max);
        pair<set<Box>::iterator,bool> ret = boxes.insert(b);
        if(ret.second)
        {
            cout<<boxes.size()<<" "<<flush;
        }
    }
    cout<<endl;
}

void Input_Generator::generate_near_boxes(Box &region, set<Box> &boxes, utype num_entries, coord_type edge, Mesh& mesh)
{
    srand((utype)time(NULL));
    iset t_ids;
    Point centroid = Point();
    Point max;
    itype rand_t_id;

    while(boxes.size()<num_entries)
    {
        rand_t_id = rand() % mesh.get_triangles_num();
        pair<iset::iterator,bool> ret1 = t_ids.insert(rand_t_id);
        if(ret1.second)
        {
            Geometry_Wrapper::get_triangle_centroid(rand_t_id,centroid,mesh);
            max.set((centroid.get_x() + edge),(centroid.get_y() + edge));
            //check that max point is into the domain
            if(!region.contains_with_all_closed_faces(max))
                continue;

            Box b = Box(centroid,max);
            pair<set<Box>::iterator,bool> ret = boxes.insert(b);
            if(ret.second)
            {
                cout<<boxes.size()<<" "<<flush;
            }
        }
    }
    cout<<endl;
}

void Input_Generator::generate_random_lines(Box &region, set<Box> &boxes, utype num_entries, coord_type edge)
{
    Point min;
    Point versor;
    Point max;

    while(boxes.size()<num_entries)
    {
        Input_Generator::generate_random_point(min,region);
        Input_Generator::generate_random_versor(versor);

        max.set((min.get_x() + versor.get_x()*edge),(min.get_y() + versor.get_y()*edge));

        //check that max point is into the domain
        if(!region.contains_with_all_closed_faces(max))
            continue;

        Box b = Box(min,max);
        pair<set<Box>::iterator,bool> ret = boxes.insert(b);
        if(ret.second)
        {
            cout<<boxes.size()<<" "<<flush;
        }
    }
    cout<<endl;
}

void Input_Generator::generate_near_lines(Box &region, set<Box> &boxes, utype num_entries, coord_type edge, Mesh& mesh)
{
    srand((utype)time(NULL));
    iset t_ids;
    Point centroid = Point();
    Point versor;
    Point max;
    itype rand_t_id;

    while(boxes.size()<num_entries)
    {
        rand_t_id = rand() % mesh.get_triangles_num();
        pair<iset::iterator,bool> ret1 = t_ids.insert(rand_t_id);
        if(ret1.second)
        {
            Geometry_Wrapper::get_triangle_centroid(rand_t_id,centroid,mesh);
            Input_Generator::generate_random_versor(versor);

            max.set((centroid.get_x() + versor.get_x()*edge),(centroid.get_y() + versor.get_y()*edge));
            //check that max point is into the domain
            if(!region.contains_with_all_closed_faces(max))
                continue;

            Box b = Box(centroid,max);
            pair<set<Box>::iterator,bool> ret = boxes.insert(b);
            if(ret.second)
            {
                cout<<boxes.size()<<" "<<flush;
            }
        }
    }
    cout<<endl;
}
