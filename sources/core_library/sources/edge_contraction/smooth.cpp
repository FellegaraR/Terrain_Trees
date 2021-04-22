/* 
 * File:   smooth.cpp
 * Author: ytsong
 * 
 * Created on April 22, 2021, 4:33 PM
 */
#include "smooth.h"

void Smooth::smooth(PRT_Tree &tree, Mesh &mesh, double lambda, int iter, cli_parameters &cli)
{

    this->iter_num = iter;
    this->lambda = lambda;
    int iterated_num =0;
    new_elevations.assign(mesh.get_vertices_num(),0);
    while(iterated_num!=iter_num){
    this->smooth_compute(tree.get_root(), mesh, tree.get_subdivision());
    this->set_zvalues(mesh);
    iterated_num++;
    }
}

void Smooth::smooth_parallel(PRT_Tree &tree, Mesh &mesh, double lambda, int iter, cli_parameters &cli)
{

    this->iter_num = iter;
    this->lambda = lambda;
    int iterated_num =0;
    new_elevations.assign(mesh.get_vertices_num(),0);
    while(iterated_num!=iter_num){
    #pragma omp parallel for
    for (int i = 0; i < tree.get_leaves_number(); i++)
    {
        Node_V *n = tree.get_leaf(i);
        if (!n->indexes_vertices())
            continue;
        this->smooth_leaf(*n,mesh);

    }
    this->set_zvalues(mesh);
    iterated_num++;
    }
}


void Smooth::smooth_compute(Node_V &n, Mesh &mesh, Spatial_Subdivision &division){

      if (n.is_leaf())
    {
        this->smooth_leaf(n,mesh);
        // 
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->smooth_compute(*n.get_son(i),mesh,division);
            }
        }
       
    }

}

void Smooth::smooth_leaf(Node_V &n, Mesh &mesh){

    if(!n.indexes_vertices())
        return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();
    itype v_range = v_end - v_start;
    leaf_VV_vec vvs(v_range,VV_vec());
    n.get_VV_vector(vvs,mesh);

    for(unsigned i=0;i<v_range;i++)
    {

        itype real_v_id=v_start+i;
        Vertex v =mesh.get_vertex(real_v_id);
        double orig_z = v.get_z();
        VV_vec &vv = vvs[i];
        coord_type zDistSum = 0;
        coord_type zSum=0;
           for(auto vid:vv)
        { 
            Vertex &v2=mesh.get_vertex(vid);
            zSum += v2.get_z();
           
        }
        coord_type new_z = lambda*(zSum/vv.size()-orig_z)+orig_z;


        new_elevations[real_v_id-1]=new_z;
        
    }
    vvs.clear();

}