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

#include "terrain_trees/node.h"
#include "basic_types/mesh.h"
#include "writer.h"

void Writer::write_node(ofstream& output, Node_T *n)
{
    if (n->is_leaf())
    {
        output << n->get_real_t_array_size();
        if (n->get_real_t_array_size() > 0)
        {
            output << endl << "  T ";
            for(RunIterator runIt = n->t_array_begin_iterator(), runEnd = n->t_array_end_iterator(); runIt != runEnd; ++runIt)
                output << *runIt << " ";
        }
    }
}


void Writer::write_node(ofstream &output, Node_V *n)
{
    if (n->is_leaf())
    {
        output << n->get_real_v_array_size() << " " << n->get_real_t_array_size();

        if (n->get_real_v_array_size() > 0)
        {
            output << endl << "  V ";
            for(RunIterator runIt = n->v_array_begin_iterator(), runEnd = n->v_array_end_iterator(); runIt != runEnd; ++runIt)
                output << *runIt << " ";
        }

        if (n->get_real_t_array_size() > 0)
        {
            output << endl << "  T ";
            for(RunIterator runIt = n->t_array_begin_iterator(), runEnd = n->t_array_end_iterator(); runIt != runEnd; ++runIt)
                output << *runIt << " ";
        }
    }
}

void Writer::write_tree_stats(IndexStatistics& indexStats)
{
    cout << indexStats.numNode << " ";
    cout << indexStats.numFullLeaf << " ";
    cout << indexStats.numEmptyLeaf << " ";
    cout << indexStats.minTreeDepth << " ";
    cout << indexStats.avgTreeDepth << " ";
    cout << indexStats.maxTreeDepth << " ";
    cout << indexStats.avg_vertices_per_leaf << " ";
    cout << indexStats.avg_completely_indexed_tri << " ";
    cout << indexStats.avg_partially_indexed_tri << " ";
    cout << indexStats.avg_overlapping_tri << " "; // different from zero only for PM criterion
    cout << indexStats.avg_leaves_for_tri << " ";
    cout << indexStats.avg_weighted_leaves_for_tri << " ";
    cout << indexStats.max_leaves_for_tri << " ";

    cout << (indexStats.numTin1Leaf*100)/(coord_type)indexStats.num_leaves_for_tri.size() << " ";
    cout << (indexStats.numTin2Leaf*100)/(coord_type)indexStats.num_leaves_for_tri.size() << " ";
    cout << (indexStats.numTin3Leaf*100)/(coord_type)indexStats.num_leaves_for_tri.size() << " ";
    cout << (indexStats.numTin4Leaf*100)/(coord_type)indexStats.num_leaves_for_tri.size() << " ";
    cout << (indexStats.numTinMoreLeaf*100)/(coord_type)indexStats.num_leaves_for_tri.size() << " ";
    cout << endl;
    // for debug only (not for tables)
    cout << "vertices_per_leaf " << indexStats.min_vertices_per_leaf << " " << indexStats.avg_vertices_per_leaf << " " << indexStats.max_vertices_per_leaf << endl;
    cout << "internal_tri_per_leaf " << indexStats.min_completely_indexed_tri << " " << indexStats.avg_completely_indexed_tri << " " << indexStats.max_completely_indexed_tri << endl;
    cout << "partial_tri_per_leaf " << indexStats.min_partially_indexed_tri << " " << indexStats.avg_partially_indexed_tri << " " << indexStats.max_partially_indexed_tri << endl;
    cout << "overlapping_tri_per_leaf " << indexStats.min_overlapping_tri << " " << indexStats.avg_overlapping_tri << " " << indexStats.max_overlapping_tri << endl;
    cout << "leaf_per_tri " << indexStats.min_leaves_for_tri << " " << indexStats.avg_leaves_for_tri << " " << indexStats.max_leaves_for_tri << " " << endl;
    cout << "chi_star " << indexStats.avg_weighted_leaves_for_tri << endl;
    cout << "t_list_length " << indexStats.t_list_length << endl;
    cout << "real_t_list_length " << indexStats.real_t_list_length << endl;
    return;
}

void Writer::write_queries_stats(int size, FullQueryStatistics& fullQueryStats, int hit_ratio)
{    
    cerr << "query_stats: ";
    cerr << fullQueryStats.minNode << " ";
    cerr << fullQueryStats.avgNode / (coord_type) size << " ";
    cerr << fullQueryStats.maxNode << " ";
    cerr << fullQueryStats.minLeaf << " ";
    cerr << fullQueryStats.avgLeaf / (coord_type) size << " ";
    cerr << fullQueryStats.maxLeaf << " ";
    cerr << fullQueryStats.mintri << " ";
    cerr << fullQueryStats.avgtri / (coord_type) size << " ";
    cerr << fullQueryStats.maxtri << " ";
    cerr << fullQueryStats.minGeometricTest << " ";
    cerr << fullQueryStats.avgGeometricTest / (coord_type) size << " ";
    cerr << fullQueryStats.maxGeometricTest << " ";

    //should be printed only for box queries
    if(fullQueryStats.maxMultipletriAccess > fullQueryStats.minMultipletriAccess)
    {
        cerr << fullQueryStats.minUniquetriAccess << " ";
        cerr << fullQueryStats.avgUniquetriAccess / (coord_type) size << " ";
        cerr << fullQueryStats.maxUniquetriAccess << " ";
        cerr << fullQueryStats.minMultipletriAccess << " ";
        cerr << fullQueryStats.avgMultipletriAccess / (coord_type) size << " ";
        cerr << fullQueryStats.maxMultipletriAccess << " ";
    }
    cerr << hit_ratio << endl;
}

void Writer::write_point_queries(set<Point>& points, string fileName)
{
    ofstream output(fileName.c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    output << points.size() << endl;
    for(set<Point>::iterator it=points.begin(); it!=points.end(); ++it)
    {
        const Point &p = *it;
        output << p.get_x() << " " << p.get_y() << endl;
    }
    output.close();
}

void Writer::write_box_queries(set<Box> &boxes, string fileName)
{
    ofstream output(fileName.c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    output << boxes.size() << endl;
    for(set<Box>::iterator it=boxes.begin(); it!=boxes.end(); ++it)
    {
        output << *it << endl;
    }
    output.close();
}

void Writer::write_mesh_VTK(string mesh_name, Mesh &mesh)
{
    stringstream stream;
    stream<<mesh_name<<".vtk";
    ofstream output(stream.str().c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    output<<"# vtk DataFile Version 2.0" << endl << endl
         << "ASCII" << endl << "DATASET UNSTRUCTURED_GRID " <<  endl << endl;

    output<< "POINTS " << mesh.get_vertices_num() << " float" << endl;

    for(itype v=1; v<=mesh.get_vertices_num(); v++)
    {
        Vertex& vert = mesh.get_vertex(v);
        output<<vert.get_x()<<" "<<vert.get_y()<<" "<<vert.get_z()<<endl;
    }

    output<<endl << "CELLS " << mesh.get_triangles_num() << " " << (mesh.get_triangles_num()*4) << endl;

    for(itype t=1; t<=mesh.get_triangles_num(); t++)
    {
        output<<"3 ";
        for(int i=0; i< mesh.get_triangle(t).vertices_num(); i++)
            output<<mesh.get_triangle(t).TV(i)-1<<" ";
        output<<endl;
    }

    output<< endl << "CELL_TYPES " << mesh.get_triangles_num() << endl;
    for (itype i = 0; i < mesh.get_triangles_num(); ++i)
        output<< "6 ";
    output<< endl;

    output<< "POINT_DATA " << mesh.get_vertices_num() << endl << endl;
    output<< "FIELD FieldData 1" << endl << endl;
    output<< "fieldvalue 1 " << mesh.get_vertices_num() << " float" << endl;

    for (itype i=1; i <=mesh.get_vertices_num(); ++i)
        output<<mesh.get_vertex(i).get_z()<< " ";
    output<<endl;

    output.close();
}

void Writer::write_mesh_curvature_VTK(string mesh_name, Mesh &mesh, string curvature_type,int c_pos/*, dvect &curvatures*/)
{
    stringstream stream;
    stream<<mesh_name<<"_"<<curvature_type<<"_curvatures_.vtk";
    ofstream output(stream.str().c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    output<<"# vtk DataFile Version 2.0" << endl << endl
         << "ASCII" << endl << "DATASET UNSTRUCTURED_GRID " <<  endl << endl;

    output<< "POINTS " << mesh.get_vertices_num() << " float" << endl;

    for(itype v=1; v<=mesh.get_vertices_num(); v++)
    {
        Vertex& vert = mesh.get_vertex(v);
        output<<vert.get_x()<<" "<<vert.get_y()<<" "<<"0"<<endl;
    }

    output<<endl << "CELLS " << mesh.get_triangles_num() << " " << (mesh.get_triangles_num()*4) << endl;

    for(itype t=1; t<=mesh.get_triangles_num(); t++)
    {
        output<<"3 ";
        for(int i=0; i< mesh.get_triangle(t).vertices_num(); i++)
            output<<mesh.get_triangle(t).TV(i)-1<<" ";
        output<<endl;
    }

    output<< endl << "CELL_TYPES " << mesh.get_triangles_num() << endl;
    for (itype i = 0; i < mesh.get_triangles_num(); ++i)
        output<< "6 ";
    output<< endl;

    output<< "POINT_DATA " << mesh.get_vertices_num() << endl << endl;
    output<< "FIELD FieldData 1" << endl << endl;
    output<< "fieldvalue 1 " << mesh.get_vertices_num() << " float" << endl;

    for(itype v=1; v<=mesh.get_vertices_num(); v++)
    {
        Vertex& vert = mesh.get_vertex(v);
        output<<vert.get_field(c_pos)<<" ";
    }
    output<<endl;

    output.close();
}



void Writer::write_mesh_roughness_VTK(string mesh_name, Mesh &mesh,int c_pos/*, dvect &curvatures*/)
{
    stringstream stream;
    stream<<mesh_name<<"_roughness_.vtk";
    ofstream output(stream.str().c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    output<<"# vtk DataFile Version 2.0" << endl << endl
         << "ASCII" << endl << "DATASET UNSTRUCTURED_GRID " <<  endl << endl;

    output<< "POINTS " << mesh.get_vertices_num() << " float" << endl;

    for(itype v=1; v<=mesh.get_vertices_num(); v++)
    {
        Vertex& vert = mesh.get_vertex(v);
        output<<vert.get_x()<<" "<<vert.get_y()<<" "<<"0"<<endl;
    }

    output<<endl << "CELLS " << mesh.get_triangles_num() << " " << (mesh.get_triangles_num()*4) << endl;

    for(itype t=1; t<=mesh.get_triangles_num(); t++)
    {
        output<<"3 ";
        for(int i=0; i< mesh.get_triangle(t).vertices_num(); i++)
            output<<mesh.get_triangle(t).TV(i)-1<<" ";
        output<<endl;
    }

    output<< endl << "CELL_TYPES " << mesh.get_triangles_num() << endl;
    for (itype i = 0; i < mesh.get_triangles_num(); ++i)
        output<< "6 ";
    output<< endl;

    output<< "POINT_DATA " << mesh.get_vertices_num() << endl << endl;
    output<< "FIELD FieldData 1" << endl << endl;
    output<< "fieldvalue 1 " << mesh.get_vertices_num() << " float" << endl;

    for(itype v=1; v<=mesh.get_vertices_num(); v++)
    {
        Vertex& vert = mesh.get_vertex(v);
        output<<vert.get_field(c_pos)<<" ";
    }
    output<<endl;

    output.close();
}


void Writer::write_mesh_gradient_VTK(string mesh_name, Mesh &mesh,int c_pos/*, dvect &curvatures*/)
{
    stringstream stream;
    stream<<mesh_name<<"_gradient_.vtk";
    ofstream output(stream.str().c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    output<<"# vtk DataFile Version 2.0" << endl << endl
         << "ASCII" << endl << "DATASET UNSTRUCTURED_GRID " <<  endl << endl;

    output<< "POINTS " << mesh.get_vertices_num() << " float" << endl;

    for(itype v=1; v<=mesh.get_vertices_num(); v++)
    {
        Vertex& vert = mesh.get_vertex(v);
        output<<vert.get_x()<<" "<<vert.get_y()<<" "<<vert.get_z()<<endl;
    }

    output<<endl << "CELLS " << mesh.get_triangles_num() << " " << (mesh.get_triangles_num()*4) << endl;

    for(itype t=1; t<=mesh.get_triangles_num(); t++)
    {
        output<<"3 ";
        for(int i=0; i< mesh.get_triangle(t).vertices_num(); i++)
            output<<mesh.get_triangle(t).TV(i)-1<<" ";
        output<<endl;
    }

    output<< endl << "CELL_TYPES " << mesh.get_triangles_num() << endl;
    for (itype i = 0; i < mesh.get_triangles_num(); ++i)
        output<< "6 ";
    output<< endl;

    output<< "POINT_DATA " << mesh.get_vertices_num() << endl << endl;
    output<< "FIELD FieldData 1" << endl << endl;
    output<< "fieldvalue 1 " << mesh.get_vertices_num() << " float" << endl;

    for(itype v=1; v<=mesh.get_vertices_num(); v++)
    {
        Vertex& vert = mesh.get_vertex(v);
        output<<"gradient_x:"<<vert.get_field(c_pos)<<" ";
    }
    output<<endl;
    
    output<< "FIELD FieldData 1" << endl << endl;
    output<< "fieldvalue 2 " << mesh.get_vertices_num() << " float" << endl;

    for(itype v=1; v<=mesh.get_vertices_num(); v++)
    {
        Vertex& vert = mesh.get_vertex(v);
        output<<"gradient_x:"<<vert.get_field(c_pos)<<" ";
    }
    output<<endl;
    output.close();
}

void Writer::write_mesh_multifield_VTK(string mesh_name, Mesh &mesh,int c_pos, string mode )
{
     stringstream stream;
    stream<<mesh_name<<"_multifield_"<<mode<<".vtk";
    ofstream output(stream.str().c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    output<<"# vtk DataFile Version 2.0" << endl << endl
         << "ASCII" << endl << "DATASET UNSTRUCTURED_GRID " <<  endl << endl;

    output<< "POINTS " << mesh.get_vertices_num() << " float" << endl;

    for(itype v=1; v<=mesh.get_vertices_num(); v++)
    {
        Vertex& vert = mesh.get_vertex(v);
        output<<vert.get_x()<<" "<<vert.get_y()<<" 0"<<endl;
    }

    output<<endl << "CELLS " << mesh.get_triangles_num() << " " << (mesh.get_triangles_num()*4) << endl;

    for(itype t=1; t<=mesh.get_triangles_num(); t++)
    {
        output<<"3 ";
        for(int i=0; i< mesh.get_triangle(t).vertices_num(); i++)
            output<<mesh.get_triangle(t).TV(i)-1<<" ";
        output<<endl;
    }

    output<< endl << "CELL_TYPES " << mesh.get_triangles_num() << endl;
    for (itype i = 0; i < mesh.get_triangles_num(); ++i)
        output<< "6 ";
    output<< endl;

    output<< "POINT_DATA " << mesh.get_vertices_num() << endl << endl;
    output<< "FIELD FieldData 1" << endl << endl;
    output<< "fieldvalue 1 " << mesh.get_vertices_num() << " float" << endl;

    for(itype v=1; v<=mesh.get_vertices_num(); v++)
    {
        Vertex& vert = mesh.get_vertex(v);
        output<<vert.get_field(c_pos)<<" ";
    }
    output<<endl;

    output.close();
}


void Writer::write_mesh(string mesh_name, string operation_type, Mesh &mesh, bool extra_fields)
{
    stringstream stream;
    stream<<mesh_name<<"_"<<operation_type<<".off";
    ofstream output(stream.str().c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    output<<"OFF"<<endl;
    output<<mesh.get_vertices_num()<<" "<<mesh.get_triangles_num()<<" 0"<<endl;
    for(itype v=1; v<=mesh.get_vertices_num(); v++)
    {
        if(!mesh.is_vertex_removed(v))
        {
            Vertex& vert = mesh.get_vertex(v);
            output<<vert.get_x()<<" "<<vert.get_y()<<" "<<vert.get_z()<<" ";
            if(extra_fields)
            {
                for(int i=1; i<vert.get_fields_num(); i++) /// NOTA: WE DO NOT DUPLICATE THE Z-COORDINATE
                    output<<vert.get_field(i)<<" ";
            }
            output<<endl;
        }
    }
    for(itype t=1; t<=mesh.get_triangles_num(); t++)
    {
        if(!mesh.is_triangle_removed(t))
        {
            Triangle &tri = mesh.get_triangle(t);
            output<<tri.vertices_num()<<" "<<tri.TV(0)-1<<" "<<tri.TV(1)-1<<" "<<tri.TV(2)-1<<endl;
        }
    }

    output.close();
}

void Writer::write_filtered_points_cloud(string mesh_name, Mesh &mesh)
{
    stringstream stream;
    stream<<mesh_name<<".xy";

    ofstream output(stream.str().c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    itype w_v_counter = 0;

    for(itype v=1; v<=mesh.get_vertices_num(); v++)
    {
        if(!mesh.is_vertex_removed(v))
        {
            Vertex& vert = mesh.get_vertex(v);
            output<<vert.get_x()<<","<<vert.get_y()<<endl;
            w_v_counter++;
        }
    }
    cout<<"[STATS] written vertices: "<<w_v_counter<<" of the initial vertices in the mesh: "<<mesh.get_vertices_num()<<endl;
    output.close();
}

void Writer::write_filtered_points_cloud_with_id(string mesh_name, Mesh &mesh)
{
    stringstream stream;
    stream<<mesh_name<<".idxy";
    ofstream output(stream.str().c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    itype w_v_counter = 0;

    for(itype v=1; v<=mesh.get_vertices_num(); v++)
    {
        if(!mesh.is_vertex_removed(v))
        {
            Vertex& vert = mesh.get_vertex(v);
            output<<v<<","<<vert.get_x()<<","<<vert.get_y()<<endl;
            w_v_counter++;
        }
    }
    cout<<"[STATS] written vertices: "<<w_v_counter<<" of the initial vertices in the mesh: "<<mesh.get_vertices_num()<<endl;
    output.close();
}

void Writer::write_multifield_points_cloud(string mesh_name, vertex_multifield &multifield, Mesh &mesh)
{
    stringstream stream;
    stream<<mesh_name<<".multi";
    ofstream output(stream.str().c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    output<<multifield.size()<<endl;
    for(vertex_multifield::iterator it=multifield.begin(); it!=multifield.end(); ++it)
    {
        Vertex& vert = mesh.get_vertex(it->first);
        output<<vert.get_x()<<" "<<vert.get_y()<<endl;

        dset &mf = it->second;
        for(dset_iter itmf=mf.begin(); itmf!=mf.end(); ++itmf)
        {
            output<<*itmf<<" ";
        }
        output<<endl;
    }
    output.close();
}

void Writer::write_critical_points(string mesh_name, map<short,set<ivect> > &critical_simplices, Mesh &mesh)
{
    stringstream stream, stream2;
    stream<<mesh_name<<"_critical_points.vtk";
    stream2<<mesh_name<<"_critical_points.txt";

    ofstream output2(stream2.str().c_str());
    output2.unsetf( std::ios::floatfield ); // floatfield not set
    output2.precision(15);

    list<Vertex > b_points;
    for(utype d=0; d<critical_simplices.size(); d++)
    {
        for(auto s : critical_simplices[d])
        {
            Vertex v = mesh.get_vertex(s[0]);

            for(auto i : s)
                output2<<mesh.get_vertex(i).get_x()<<" "<<mesh.get_vertex(i).get_y()<<" "<<mesh.get_vertex(i).get_z()<<" - ";
            output2<<endl;

            for(itype i=1; i<s.size(); i++)
            {
                Vertex &v2 = mesh.get_vertex(s[i]);
                for(int i=0; i<=v.get_dimension(); i++)
                {
                    v.set_c(i,((v.get_c(i)+v2.get_c(i)) / 2.0));
                }
            }
            b_points.push_back(v);
        }
    }
    output2.close();

    ofstream output(stream.str().c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    output<< "# vtk DataFile Version 2.0"<<endl<<endl;
    output<< "ASCII "<<endl;
    output<< "DATASET UNSTRUCTURED_GRID"<<endl<<endl;;
    output<< "POINTS "<<b_points.size()<<" float"<<endl;

    for(auto v : b_points)
        output<<v.get_x()<<" "<<v.get_y()<<" "<<v.get_z()<<endl;
    output<<endl<<endl;

    output<< "POINT_DATA "<<b_points.size()<<endl;
    output<< "FIELD FieldData 1"<<endl;
    output<< "cp_type 1 "<<b_points.size()<<" float"<<endl;

    for(utype d=0; d<critical_simplices.size(); d++)
    {
        for(auto s : critical_simplices[d])
            output<<s.size()-1<<" ";
    }
    output<<endl;
    output.close();
}

void Writer::write_field_csv(string mesh_name, Mesh& mesh){

   stringstream stream;
    stream<<mesh_name<<"_field.csv";
    ofstream output(stream.str().c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

   

    output<< "x,y,z,class,R,G,B,roughness,multifield"<< endl;
    
    for(itype v=1; v<=mesh.get_vertices_num(); v++)
    {
        Vertex& vert = mesh.get_vertex(v);
        output<<vert.get_x()<<","<<vert.get_y()<<","<<vert.get_z();
        for(itype i=1;i<=6;i++)
        {
            output<<","<<vert.get_field(i);
        }
        output<<endl;
       
    }

    output.close();
}
