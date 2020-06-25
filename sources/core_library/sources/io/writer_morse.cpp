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

#include "writer_morse.h"

void Writer_Morse::write_asc1cells_VTK(string mesh_name, string operation_type, itype vertices_per_leaf, simplices_map &triangles, Mesh &mesh,
                                  ivect &original_triangle_indices, ivect &original_vertex_indices, dvect &original_vertex_fields, bool revert_to_original_field)
{
    //init recupero i vertici unici e li ordino
    ivect new_vertex_index(mesh.get_vertices_num(), -1);
    ivect vertici_ordinati;

    itype vertex_number=0;
    for(simplices_map::iterator it = triangles.begin(); it!=triangles.end(); it++)
    {
        for(int i=0; i<3; i++)
        {
            if(new_vertex_index[(it->first)[i]-1] == -1)
            {
                new_vertex_index[(it->first)[i]-1] = vertex_number++;
                vertici_ordinati.push_back((it->first)[i]);
            }
        }
    }
    //fine init

    stringstream stream;
    stream<<mesh_name<<"_kv_"<<vertices_per_leaf<<"_"<<operation_type<<".vtk";
    ofstream output(stream.str().c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    output<<"# vtk DataFile Version 2.0" << endl << endl
         << "ASCII" << endl << "DATASET UNSTRUCTURED_GRID " <<  endl << endl;

    output<< "POINTS " << vertex_number << " float" << endl;

    for(itype i=0; i<vertici_ordinati.size(); i++)
    {
        Vertex& vert = mesh.get_vertex(vertici_ordinati.at(i));
        for(int i=0; i<vert.get_dimension(); i++)
            output<<vert.get_c(i)<<" ";
        if(vert.get_dimension()==2 && revert_to_original_field)
            output<<original_vertex_fields.at(original_vertex_indices.at(vertici_ordinati.at(i)-1));
        else
            output<<vert.get_z();
        output<<endl;
    }
    output<<endl;
    output<<endl;

    output<<endl << "CELLS " << triangles.size() << " " << (triangles.size()*4) << endl;

    for(simplices_map::iterator it = triangles.begin(); it != triangles.end(); it++)
    {
        output<<"3 "<<new_vertex_index[(it->first)[0]-1]<<" "<<new_vertex_index[(it->first)[1]-1]<<" "<<new_vertex_index[(it->first)[2]-1]<<endl;
    }
    output<<endl;

    output<< endl << "CELL_TYPES " << triangles.size() << endl;
    for (itype i = 0; i < triangles.size(); ++i)
        output<< "5 ";
    output<< endl;
    output<< endl;

    output<< "POINT_DATA " << vertex_number << endl << endl;
    output<< "FIELD FieldData 1" << endl << endl;
    output<< "original_field 1 " << vertex_number << " float" << endl;

    for (itype v=1; v <=mesh.get_vertices_num(); ++v)
    {
        if(new_vertex_index[v-1] != -1)
        {
            if(revert_to_original_field)
                output<<original_vertex_fields.at(original_vertex_indices.at(v-1))<<" ";
            else
                output<<mesh.get_vertex(v).get_z()<< " ";
        }
    }

    output<<endl;
    output<<endl;

    output<< "CELL_DATA " << triangles.size() << endl;
    output<< "FIELD FieldData 1" << endl << endl;
    if(operation_type == "asc1cells")
        output<< "ascending_1_cells 1 " << triangles.size() << " int" << endl;
    else
        output<< "descending_2_cells 1 " << triangles.size() << " int" << endl;

    for(simplices_map::iterator it = triangles.begin(); it != triangles.end(); it++)
    {
        output << original_triangle_indices.at(it->second-1) << " ";
    }
    output<<endl;

    output.close();
}


void Writer_Morse::write_desc2cells_VTK(string mesh_name, string operation_type, itype vertices_per_leaf, ivect &segmentation, Mesh &mesh,
                                  ivect &original_triangle_indices, ivect &original_vertex_indices, dvect &original_vertex_fields, bool revert_to_original_field)
{
    stringstream stream;
    stream<<mesh_name<<"_kv_"<<vertices_per_leaf<<"_"<<operation_type<<".vtk";
    ofstream output(stream.str().c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    output<<"# vtk DataFile Version 2.0" << endl << endl
         << "ASCII" << endl << "DATASET UNSTRUCTURED_GRID " <<  endl << endl;

    output<< "POINTS " << mesh.get_vertices_num() << " float" << endl;

    for(itype v=1; v<=mesh.get_vertices_num(); v++)
    {
        Vertex& vert = mesh.get_vertex(v);
        for(int i=0; i<vert.get_dimension(); i++)
            output<<vert.get_c(i)<<" ";
        if(vert.get_dimension()==2 && revert_to_original_field)
            output<<original_vertex_fields.at(original_vertex_indices.at(v-1));
        else
            output<<vert.get_z();
        output<<endl;
    }

    output<<endl << "CELLS " << mesh.get_triangles_num() << " " << (mesh.get_triangles_num()*4) << endl;

    for(itype t=1; t<=mesh.get_triangles_num(); t++)
    {
        output<<"3 ";
        for(int i=0; i< mesh.get_triangle(t).vertices_num(); i++)
            output<<abs(mesh.get_triangle(t).TV(i))-1<<" ";
        output<<endl;
    }

    output<< endl << "CELL_TYPES " << mesh.get_triangles_num() << endl;
    for (itype i = 0; i < mesh.get_triangles_num(); ++i)
        output<< "5 ";
    output<< endl;

    output<< "POINT_DATA " << mesh.get_vertices_num() << endl << endl;
    output<< "FIELD FieldData 1" << endl << endl;
    output<< "original_field 1 " << mesh.get_vertices_num() << " float" << endl;

    for (itype v=1; v <=mesh.get_vertices_num(); ++v)
    {
        if(revert_to_original_field)
            output<<original_vertex_fields.at(original_vertex_indices.at(v-1))<<" ";
        else
            output<<mesh.get_vertex(v).get_z()<< " ";
    }
    output<<endl;
    output<<endl;

    output<< "CELL_DATA " << mesh.get_triangles_num() << endl;
    output<< "FIELD FieldData 1" << endl << endl;
    output<< "descending_2_cells 1 " << mesh.get_triangles_num() << " float" << endl;

    for(itype i=0; i<mesh.get_triangles_num(); i++)
    {
        if(segmentation[i] != -1)
            output<< original_triangle_indices.at(segmentation[i]-1) << " ";
        else
            output << segmentation[i] << " ";
    }

    output<<endl;

    output.close();
}

void Writer_Morse::write_desc1cells_VTK(string mesh_name, string operation_type, itype vertices_per_leaf, simplices_map &edges, Mesh &mesh,
                                  ivect &original_vertex_indices, dvect &original_vertex_fields, bool revert_to_original_field)
{
    itype edge_number = edges.size();

    ivect new_vertex_index = ivect(mesh.get_vertices_num(), -1);
    ivect vertici_ordinati;

    itype vertex_number=0;
    for(simplices_map::iterator it = edges.begin(); it!=edges.end(); ++it)
    {
        for(uint i=0; i<(it->first).size(); i++)
        {
            if(new_vertex_index[(it->first).at(i)-1] == -1)
            {
                new_vertex_index[(it->first).at(i)-1] = vertex_number++;
                vertici_ordinati.push_back((it->first).at(i));
            }
        }
    }

    stringstream stream;
    stream<<mesh_name<<"_kv_"<<vertices_per_leaf<<"_"<<operation_type<<".vtk";
    ofstream output(stream.str().c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    output<<"# vtk DataFile Version 2.0" << endl << endl
         << "ASCII" << endl << "DATASET UNSTRUCTURED_GRID " <<  endl << endl;

    output<< "POINTS " << vertex_number << " float" << endl;

    for(itype i=0; i<vertici_ordinati.size(); i++)
    {
        Vertex& vert = mesh.get_vertex(vertici_ordinati.at(i));
        for(int i=0; i<vert.get_dimension(); i++)
            output<<vert.get_c(i)<<" ";
        if(vert.get_dimension()==2 && revert_to_original_field)
            output<<original_vertex_fields.at(original_vertex_indices.at(vertici_ordinati.at(i)-1));
        else
            output<<vert.get_z();
        output<<endl;
    }
    output<<endl;
    output<<endl;

    output<<endl << "CELLS " << edge_number << " " << (edge_number*3) << endl;

    for(simplices_map::iterator it = edges.begin(); it != edges.end(); it++)
    {
        output<<"2 "<<new_vertex_index[it->first[0]-1]<<" "<<new_vertex_index[it->first[1]-1]<<endl;

        if(new_vertex_index[it->first[0]-1] == -1 || new_vertex_index[it->first[1]-1] == -1)
        {
            cout<<"2 "<<new_vertex_index[it->first[0]-1]<<" "<<new_vertex_index[it->first[1]-1]<<endl;
            int a; cin>>a;
        }
    }
    output<<endl;

    output<< endl << "CELL_TYPES " << edge_number << endl;
    for (itype i = 0; i < edges.size(); ++i)
        output<< "3 ";
    output<< endl;
    output<< endl;

    output<< "POINT_DATA " << vertex_number << endl << endl;
    output<< "FIELD FieldData 1" << endl << endl;
    output<< "original_field 1 " << vertex_number << " float" << endl;

    for(itype i=0; i<vertici_ordinati.size(); i++)
    {
        if(revert_to_original_field)
            output<<original_vertex_fields.at(original_vertex_indices.at(vertici_ordinati.at(i)-1))<<" ";
        else
            output<<mesh.get_vertex(vertici_ordinati.at(i)).get_z()<< " ";
    }
    output<< endl;
    output<< endl;

    output<< "CELL_DATA " << edges.size() << endl;
    output<< "FIELD FieldData 1" << endl << endl;
    output<< "descending_1_cells 1 " << edges.size() << " int" << endl;

    for(simplices_map::iterator it = edges.begin(); it != edges.end(); it++)
    {
        output << original_vertex_indices.at(it->second-1)-1 << " ";
    }
    output<<endl;

    output.close();
}

void Writer_Morse::write_asc2cells_VTK(string mesh_name, string operation_type, itype vertices_per_leaf, ivect &segmentation, Mesh &mesh,
                                 ivect &original_vertex_indices, dvect &original_vertex_fields, bool revert_to_original_field)
{
    stringstream stream;
    stream<<mesh_name<<"_kv_"<<vertices_per_leaf<<"_"<<operation_type<<".vtk";
    ofstream output(stream.str().c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    output<<"# vtk DataFile Version 2.0" << endl << endl
         << "ASCII" << endl << "DATASET UNSTRUCTURED_GRID " <<  endl << endl;

    output<< "POINTS " << mesh.get_vertices_num() << " float" << endl;

    for(itype v=1; v<=mesh.get_vertices_num(); v++)
    {
        Vertex& vert = mesh.get_vertex(v);
        for(int i=0; i<vert.get_dimension(); i++)
            output<<vert.get_c(i)<<" ";
        if(vert.get_dimension()==2 && revert_to_original_field)
            output<<original_vertex_fields.at(original_vertex_indices.at(v-1));
        else
            output<<vert.get_z();
        output<<endl;
    }

    output<<endl << "CELLS " << mesh.get_triangles_num() << " " << (mesh.get_triangles_num()*4) << endl;

    for(itype t=1; t<=mesh.get_triangles_num(); t++)
    {
        output<<"3 ";
        for(int i=0; i< mesh.get_triangle(t).vertices_num(); i++)
            output<<abs(mesh.get_triangle(t).TV(i))-1<<" ";
        output<<endl;
    }

    output<< endl << "CELL_TYPES " << mesh.get_triangles_num() << endl;
    for (itype i = 0; i < mesh.get_triangles_num(); ++i)
        output<< "5 ";
    output<< endl;

    output<< "POINT_DATA " << mesh.get_vertices_num() << endl << endl;
    output<< "FIELD FieldData 2" << endl << endl;
    output<< "original_field 1 " << mesh.get_vertices_num() << " float" << endl;

    for (itype v=1; v <=mesh.get_vertices_num(); ++v)
    {
        if(revert_to_original_field)
            output<<original_vertex_fields.at(original_vertex_indices.at(v-1))<<" ";
        else
            output<<mesh.get_vertex(v).get_z()<< " ";
    }
    output<<endl;
    output<<endl;

    output<< "ascending_2_cells 1 " << mesh.get_vertices_num() << " float" << endl;

    for(itype i=0; i<mesh.get_vertices_num(); i++)
    {
        if(segmentation[i] != -1)
            output<< original_vertex_indices.at(segmentation[i]) << " ";
        else
            output << segmentation[i] << " ";
    }

    output<<endl;

    output.close();
}

void Writer_Morse::write_incidence_graph_VTK(string mesh_name, string operation_type, itype vertices_per_leaf, IG &forman_ig, Mesh &mesh,
                                       ivect &original_vertex_indices, dvect &original_vertex_fields, bool revert_to_original_field)
{
    stringstream stream;
    stream<<mesh_name<<"_kv_"<<vertices_per_leaf<<"_"<<operation_type<<".vtk";
    ofstream output(stream.str().c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    ivect new_vertex_index = ivect(mesh.get_vertices_num(), -1);
    vector<bool> connected = vector<bool>(mesh.get_vertices_num(), false);
    set<pair<itype,itype> > arcs = set<pair<itype,itype> >();

    itype vertex_number = 0;
    itype edge_number =0;

    ivect critici = ivect(mesh.get_vertices_num(), -1);

    map<itype, nNode*>& minima = forman_ig.getMinima();
    for(map<itype, nNode*>::iterator it=minima.begin(); it!=minima.end(); ++it)
    {
        nNode* node = it->second;
        critici[node->get_critical_index()-1]=0;
    }

    map<itype, nNode*>& maxima = forman_ig.getMaxima();
    for(map<itype, nNode*>::iterator it=maxima.begin(); it!=maxima.end(); ++it)
    {
        nNode* node = it->second;
        itype max_v = mesh.get_max_elevation_vertex(mesh.get_triangle(node->get_critical_index()));
        critici[max_v-1]=2;
    }

    map<pair<itype,itype>, iNode*>& saddle = forman_ig.getSaddles();
    for(map<pair<itype,itype>, iNode*>::iterator it=saddle.begin(); it!=saddle.end(); ++it)
    {
        iNode* node = it->second;
        itype sad_vertex = node->get_critical_index();
        critici[sad_vertex-1]=1;

        set<Arc*> &arcs_up = node->getArcs(true);
        for(set<Arc*>::const_iterator it=arcs_up.begin(); it!=arcs_up.end(); ++it)
        {
            arcs.insert(pair<itype,itype>(((*it)->getNode_i())->get_critical_index() ,sad_vertex));
        }

        set<Arc*> &arcs_down = node->getArcs(false);
        for(set<Arc*>::const_iterator it=arcs_down.begin(); it!=arcs_down.end(); ++it)
        {
            itype max_v = mesh.get_max_elevation_vertex(mesh.get_triangle(((*it)->getNode_j())->get_critical_index()));
            arcs.insert(pair<itype,itype>(sad_vertex, max_v));
        }
    }

    for(set<pair<itype,itype> >::iterator it = arcs.begin(); it!=arcs.end(); it++){
        connected[it->first-1] = true;
        connected[it->second-1] = true;
    }

    for(itype i=1; i<=mesh.get_vertices_num(); i++)
    {
        if(critici[i-1] != -1 && connected[i-1])
        {
            new_vertex_index[i-1] = vertex_number;
            vertex_number++;
        }
    }

    edge_number = arcs.size();

    output<<"# vtk DataFile Version 2.0" << endl << endl
         << "ASCII" << endl << "DATASET UNSTRUCTURED_GRID " <<  endl << endl;
    output<< "POINTS " << vertex_number << " float" << endl;

    for (itype v=1; v <=mesh.get_vertices_num(); ++v)
    {
        if(new_vertex_index[v-1] != -1)
        {
            Vertex& vert = mesh.get_vertex(v);
            for(int i=0; i<vert.get_dimension(); i++)
                output<<vert.get_c(i)<<" ";
            if(vert.get_dimension()==2 && revert_to_original_field)
                output<<original_vertex_fields.at(original_vertex_indices.at(v-1));
            else
                output<<vert.get_z();
            output<<endl;
        }
    }
    output<<endl;

    output<<endl << "CELLS " << edge_number << " " << (edge_number*3) << endl;

    for(set<pair<itype,itype> >::iterator it = arcs.begin(); it != arcs.end(); it++)
    {
        output<<"2 "<<new_vertex_index[it->first-1]<<" "<<new_vertex_index[it->second-1]<<endl;
    }

    output<< endl << "CELL_TYPES " << edge_number << endl;
    for (utype i = 0; i < arcs.size(); ++i)
        output<< "3 ";
    output<< endl;
    output<< endl;

    output<< "POINT_DATA " << vertex_number << /*endl <<*/ endl;
    output<< "FIELD FieldData 2" << /*endl <<*/ endl;
    output<< "original_field 1 " << vertex_number << " float" << endl;

    for (itype v=1; v <=mesh.get_vertices_num(); ++v)
    {
        if(new_vertex_index[v-1] != -1 && connected[v-1])
        {
            if(revert_to_original_field)
                output<<original_vertex_fields.at(original_vertex_indices.at(v-1))<<" ";
            else
                output<<mesh.get_vertex(v).get_z()<< " ";
        }
    }
    output<< endl;
    output<< endl;

    output<< "critical_point 1 " << vertex_number << " int" << endl;

    for (utype i=0; i <critici.size(); ++i)
        if(critici[i] != -1 && connected[i])
            output<< critici[i] << " ";

    output<<endl;
    output<<endl;

    output.close();
}

void Writer_Morse::write_critical_clusters(string mesh_name, forman_aux_structures::critical_clusters &cc, Mesh &mesh)
{
    stringstream stream;
    stream<<mesh_name<<"_critical_clusters.vtk";
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
    output<< "FIELD FieldData 2" << endl << endl;
    output<< "fieldvalue 1 " << mesh.get_vertices_num() << " float" << endl;

    for (itype i=1; i <=mesh.get_vertices_num(); ++i)
        output<<mesh.get_vertex(i).get_z()<< " ";
    output<<endl;

    output<< "minClusters 1 " << mesh.get_vertices_num() << " float" << endl;
    for (itype i=1; i <=mesh.get_vertices_num(); ++i)
        output<<cc.get_v_label(i)<< " ";
    output<<endl;

    output<< endl << "CELL_DATA " << mesh.get_triangles_num() << endl;
    output<< "FIELD FieldData 1" << endl << endl;
    output<< "maxClusters 1 " << mesh.get_triangles_num() << " float" << endl;
    for (itype i = 1; i <= mesh.get_triangles_num(); ++i)
        output<<cc.get_t_label(i)<<" ";
    output<<endl;

    output.close();
}
