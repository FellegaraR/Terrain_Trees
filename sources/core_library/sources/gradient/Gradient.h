/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Gradient.h
 * Author: ytsong
 *
 * Created on January 10, 2019, 10:35 PM
 */

#ifndef GRADIENT_H
#define GRADIENT_H
#include "terrain_trees/node_v.h"
#include "terrain_trees/node_t.h"
#include "terrain_trees/spatial_subdivision.h"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <complex>
#include <limits>
#include "utilities/timer.h"

using namespace Eigen;

class Gradient {
public:
    Gradient();
    Gradient(string mode){   
        this->mode=mode;
        if(mode=="All"){
        fields.push_back(0);
        fields.push_back(1);
        fields.push_back(2);
        fields.push_back(3);
        fields.push_back(4);
      
        }
        else if (mode=="rRGB")
        {
        
        fields.push_back(4);
        fields.push_back(1);
        fields.push_back(2);
        fields.push_back(3);
        }
        else if(mode=="rGray")
        {
            
            fields.push_back(4);
            fields.push_back(5);
        }
        else if(mode=="RGB")
        {
            
            fields.push_back(1);
            fields.push_back(2);
            fields.push_back(3);
           
        }
        else if(mode=="RGB")
        {
            
            fields.push_back(1);
            fields.push_back(2);
            fields.push_back(3);
           
        }
    };
    virtual ~Gradient();
    

    
    void VT_relation(Node_V& n, Mesh& mesh, Spatial_Subdivision &division);
    void VT_relation(Node_T &n, Box &n_dom, int level, Mesh &mesh, Spatial_Subdivision &division);
    
    void VV_relation(Node_V& n, Mesh& mesh, Spatial_Subdivision &division);
    void VV_relation(Node_T &n, Box &n_dom, int level, Mesh &mesh, Spatial_Subdivision &division);
    
    void multi_field(Node_V& n, Mesh& mesh, Spatial_Subdivision &division);
    void multi_field(Node_T &n, Box &n_dom, int level, Mesh &mesh, Spatial_Subdivision &division);
 
    
    
    inline void compute_field_stats(Mesh &mesh,int factor)
    {
    // Using mode to get the field_pos, store them in to a vector
      //  ivect fields=this->fields;

        for(auto f:fields){
        coord_type min=INFINITY, max=-INFINITY;
        
        for(itype v=1; v<=mesh.get_vertices_num(); v++)
        {
            coord_type c = mesh.get_vertex(v).get_field(f);

            if(c < min)
                min = c;
            if(c > max)
                max = c;
          //  avg += c;
        }
        coord_type range=(max-min)/factor;
        this->ranges.push_back(range);
        
        }
        
      //  }
      //  else
       //     this->range=this->range/factor;
    }
    
    inline void print_multifield_stats(Mesh &mesh, int c_pos)
    {
        coord_type min=INFINITY, max=-INFINITY, avg=0;
        for(itype v=1; v<=mesh.get_vertices_num(); v++)
        {
            coord_type c = mesh.get_vertex(v).get_field(c_pos);

            if(c < min)
                min = c;
            if(c > max)
                max = c;
            avg += c;
        }
        cerr<<"[STATS] multifield min: "<<min<<" avg: "<<avg/(coord_type)mesh.get_vertices_num()<<" max: "<<max<<endl;
    }

private:

    dvect PCE_compute(Triangle &t, Mesh &mesh, int field_index);


    FG AGS_compute(itype real_v_id,VT &vt, Mesh& mesh, int field_index);

    coord_type cross_2d(dvect i,dvect j);
    coord_type dot_2d(dvect i,dvect j);

    void multi_field_leaf(Node_V& n, Mesh& mesh);
    void multi_field_leaf(Node_T& n, Box &dom, Mesh& mesh);

    void VV_relation_leaf(Node_V& n, Mesh& mesh);
    void VV_relation_leaf(Node_T& n, Box &dom, Mesh& mesh);
    
    void VT_relation_leaf(Node_V& n, Mesh& mesh);
    void VT_relation_leaf(Node_T& n, Box &dom, Mesh& mesh);
    
    
    coord_type block_time;


dvect ranges;
ivect fields;
string mode;
leaf_VT vts;
leaf_VV vvs;
};

#endif /* GRADIENT_H */

