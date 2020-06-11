/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Gradient.cpp
 * Author: ytsong
 * 
 * Created on January 10, 2019, 10:35 PM
 */

#include "Gradient.h"

Gradient::Gradient() {

}


Gradient::~Gradient() {
}


    
    dvect Gradient::PCE_compute(Triangle& t, Mesh& mesh, int field_index){
    
            coord_type area;
          
            Vertex &vi = mesh.get_vertex(t.TV(0));
            Vertex &vj = mesh.get_vertex(t.TV(1));
            Vertex &vk = mesh.get_vertex(t.TV(2));

            dvect ki={vi.get_x()-vk.get_x(),vi.get_y()-vk.get_y()};
            dvect ij={vj.get_x()-vi.get_x(),vj.get_y()-vi.get_y()};
                // compute the vectors of the edges after rotate by 90 degrees
            dvect ki_vert= { vk.get_y()-vi.get_y() , vi.get_x()-vk.get_x()};
            dvect ij_vert= {vi.get_y()-vj.get_y() , vj.get_x()-vi.get_x()};
            // compute the area of the triangle 
            area=abs(0.5*(ki[0]*ij[1]-ki[1]*ij[0]));
            // using the z coordinate as function value
            int field_pos=fields[field_index];
           /* if(mode=="rGray"&&field_index==0){
                coord_type vj_grayscale=0.2126*vj.get_field(field_pos)+0.7152*vj.get_field(field_pos+1)+0.0722*vj.get_field(field_pos+2);
                coord_type vi_grayscale=0.2126*vi.get_field(field_pos)+0.7152*vi.get_field(field_pos+1)+0.0722*vi.get_field(field_pos+2);
                coord_type vk_grayscale=0.2126*vk.get_field(field_pos)+0.7152*vk.get_field(field_pos+1)+0.0722*vk.get_field(field_pos+2);
                
                dvect area_multiply_ki={(vj_grayscale-vi_grayscale)*ki_vert[0]/(2*area),(vj_grayscale-vi_grayscale)*ki_vert[1]/(2*area)};
                dvect area_multiply_ij={(vk_grayscale-vi_grayscale)*ij_vert[0]/(2*area),(vk_grayscale-vi_grayscale)*ij_vert[1]/(2*area)};
                coord_type range=65280/100;
                dvect g={(area_multiply_ki[0]+area_multiply_ij[0])/range,(area_multiply_ki[1]+area_multiply_ij[1])/range};
                return g;
            }*/
         //   else{
            dvect area_multiply_ki={(vj.get_field(field_pos)-vi.get_field(field_pos))*ki_vert[0]/(2*area),(vj.get_field(field_pos)-vi.get_field(field_pos))*ki_vert[1]/(2*area)};
            dvect area_multiply_ij={(vk.get_field(field_pos)-vi.get_field(field_pos))*ij_vert[0]/(2*area),(vk.get_field(field_pos)-vi.get_field(field_pos))*ij_vert[1]/(2*area)};
            coord_type range=ranges[field_index];
            dvect g={(area_multiply_ki[0]+area_multiply_ij[0])/range,(area_multiply_ki[1]+area_multiply_ij[1])/range};
           // cout<<"the gradient of this triangle is "<<g[0]<<","<<g[1]<<endl;
            return g;//}
    }
    

    FG Gradient::AGS_compute(itype real_v_id,VT &vt, Mesh& mesh, int field_index){
        
        Vertex &vp=mesh.get_vertex(real_v_id);// vp is always the studied vertex
        coord_type sum_area=0.0;
        dvect sum_gradient={0,0};
        dvect gradient_v={0,0};
        for(auto tid:vt)
        {//compute the area of each triangle
            dvect gradient_t;
        Triangle &t=mesh.get_triangle(tid);
        Vertex &vq=(t.TV(0)!=real_v_id)?mesh.get_vertex(t.TV(0)):mesh.get_vertex(t.TV(1));
        Vertex &vr=((t.TV(1)!=real_v_id)&&(t.TV(0)!=real_v_id))?mesh.get_vertex(t.TV(1)):mesh.get_vertex(t.TV(2));

               gradient_t=this->PCE_compute(t,mesh,field_index);

            dvect qp={vq.get_x()-vp.get_x(),vq.get_y()-vp.get_y()};
            dvect qr={vq.get_x()-vr.get_x(),vq.get_y()-vr.get_y()};
            dvect rp={vr.get_x()-vp.get_x(),vr.get_y()-vp.get_y()};
            dvect rq={-qr[0],-qr[1]};
          //  dvect pr={-rp[0],-rp[1]};
            // compute the area of the triangle 
            coord_type Voronoi_area=0;
            if(dot_2d(qp,rp)<=0)
            { Voronoi_area=abs(0.25*(qp[0]*rp[1]-qp[1]*rp[0]));
            }
            else if((dot_2d(rp,rq)<0)||(dot_2d(qp,qr)<0))
            {Voronoi_area=abs(0.125*(qp[0]*rp[1]-qp[1]*rp[0]));
          //  cout<<"obtuse at Q: "<<dot_2d(qp,qr)<<"obtuse at R: "<<dot_2d(rp,rq)<<endl;
            }
            else            
            {Voronoi_area=0.125*(abs((rp[0]*rp[0]+rp[1]*rp[1])*dot_2d(qp,qr)/(cross_2d(qp,qr)))+abs((qp[0]*qp[0]+qp[1]*qp[1])*dot_2d(rp,rq)/cross_2d(rp,rq)));
           // cout<<"no obtuse"<<endl;
            }
         //  cout<<"Voronoi Area: "<<Voronoi_area<<endl;
            sum_area+=Voronoi_area;
            
            sum_gradient[0]+=gradient_t[0]*Voronoi_area;
            sum_gradient[1]+=gradient_t[1]*Voronoi_area;
        }
        gradient_v[0]=sum_gradient[0]/sum_area;
        gradient_v[1]=sum_gradient[1]/sum_area;

        
        FG gradient;
        gradient.push_back(gradient_v[0]);
        gradient.push_back(gradient_v[1]);
        return gradient;
    }
    
    
    
    coord_type Gradient::dot_2d(dvect i,dvect j){
        return (i[0]*j[0]+i[1]*j[1]);
    }
    
    coord_type Gradient::cross_2d(dvect i,dvect j){
        return (i[0]*j[1]-i[1]*j[0]);
    }
    void Gradient::multi_field(Node_V& n, Mesh& mesh, Spatial_Subdivision &division){
        
        /// if there are no vertices in the leaf we have nothing to do..
      if (n.is_leaf())
    {
        this->multi_field_leaf(n,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->multi_field(*n.get_son(i),mesh,division);
            }
        }
    }
        
    }
    
    
    
    
    
    void Gradient::multi_field(Node_T &n, Box &n_dom, int level, Mesh &mesh, Spatial_Subdivision &division){
        
            if (n.is_leaf())
    {
        this->multi_field_leaf(n,n_dom,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(n_dom,level,i);
            int son_level = level +1;
            if(n.get_son(i)!=NULL)
            {
                this->multi_field(*n.get_son(i),son_dom,son_level,mesh,division);
            }
        }
    }
    }
    
    
    void Gradient::multi_field_leaf(Node_V& n, Mesh& mesh){
         if(!n.indexes_vertices())
        return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();
    itype v_range = v_end - v_start;
    
    leaf_VT vts(v_range,VT());
                           Timer time;
        time.start();
    n.get_VT(vts,mesh);
      time.stop();
               block_time+=time.get_elapsed_time();
        for(unsigned i=0;i<v_range;i++)
    {
       itype real_v_id=v_start+i;
        VT vt=vts[i];
      
        vect_FG FG_matrix;
         

        for(int i=0;i<fields.size();i++){

            FG_matrix.push_back(AGS_compute(real_v_id,vt,  mesh, i));  

        }
            
        itype field_num=FG_matrix.size();
        MatrixX2d gradientMatrix(field_num,2);
        for(int i=0;i<field_num;i++)
            for(int j=0;j<2;j++)
            {
                gradientMatrix(i,j)=FG_matrix[i][j];

            }        

        MatrixXd input=gradientMatrix.transpose()*gradientMatrix;
        EigenSolver<MatrixXd> solver;
        solver.compute(input,false);   
        double max = numeric_limits<double>::min();
	for(int i=0; i<solver.eigenvalues().rows(); i++)
	{
		complex<double> c = solver.eigenvalues().coeff(i,0);
        	// cout << abs(c) << endl;
		if(max < abs(c))
			max = abs(c);
	}
 //       cout<<"the result is:"<<sqrt(max)<<endl;
        coord_type multifield=0;
        multifield=sqrt(max);
        mesh.get_vertex(real_v_id).add_field(multifield);
        }
    
    
    
    
    }
    
    void Gradient::multi_field_leaf(Node_T& n, Box &dom, Mesh& mesh){
        
    itype v_start;
    itype v_end;

    n.get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;
    itype v_range=v_end-v_start;
     
    leaf_VT vts(v_end-v_start,VT());
    n.get_VT(vts,v_start,v_end,mesh);
      for(unsigned i=0;i<v_range;i++)
    {
       itype real_v_id=v_start+i;

        VT &vt=vts[i];
        //to-do:using the mode to decide what gradient to be included.
        // here just we just use All mode.        
        
        vect_FG FG_matrix;
        for(int i=0;i<fields.size();i++){
            //cout<<"Add one field"<<endl;
            FG_matrix.push_back(AGS_compute(real_v_id,vt, mesh, i));
        
        }
  
        itype field_num=FG_matrix.size();
     //   cout<<"field num:"<<field_num<<endl;
        MatrixX2d gradientMatrix(field_num,2);
     
        for(int i=0;i<field_num;i++)
            for(int j=0;j<2;j++)
            {
                gradientMatrix(i,j)=FG_matrix[i][j];
            }        
        
        MatrixXd input=gradientMatrix.transpose()*gradientMatrix;
        EigenSolver<MatrixXd> solver;
        solver.compute(input,false);
   //     cout<<"here is the eigenvalue:"<<solver.eigenvalues()<<endl;
   
        double max = numeric_limits<double>::min();
	for(int i=0; i<solver.eigenvalues().rows(); i++)
	{
		complex<double> c = solver.eigenvalues().coeff(i,0);
	//	 cout << abs(c) << endl;
		if(max < abs(c))
			max = abs(c);
	}
  //      cout<<"the result is:"<<sqrt(max)<<endl;
        coord_type multifield=0;
        multifield=sqrt(max);
        mesh.get_vertex(real_v_id).add_field(multifield);
        }
    
    
    }
    
    
    void Gradient::VT_relation(Node_V& n, Mesh& mesh, Spatial_Subdivision& division){
    
                /// if there are no vertices in the leaf we have nothing to do..
      if (n.is_leaf())
    {
        this->VT_relation_leaf(n,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->VT_relation(*n.get_son(i),mesh,division);
            }
        }
    }
    
    }
    
    void Gradient::VT_relation(Node_T &n, Box &n_dom, int level, Mesh &mesh, Spatial_Subdivision &division){
        
            if (n.is_leaf())
    {
        this->VT_relation_leaf(n,n_dom,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(n_dom,level,i);
            int son_level = level +1;
            if(n.get_son(i)!=NULL)
            {
                this->VT_relation(*n.get_son(i),son_dom,son_level,mesh,division);
            }
        }
    }
    }
    
    void Gradient::VT_relation_leaf(Node_V& n, Mesh& mesh){
    
    //    if(!n.indexes_vertices())
    //     return;

    // itype v_start = n.get_v_start();
    // itype v_end = n.get_v_end();
    // itype v_range = v_end - v_start;
    
    leaf_VT vts;//(v_range,VT());

    n.get_VT(vts,mesh);
   
    //     for(unsigned i=0;i<v_range;i++)
    // {
    //    itype real_v_id=v_start+i;
    //     VT &vt=vts[i];
      
    //  }  
    
    }
    
    void Gradient::VT_relation_leaf(Node_T& n, Box &dom, Mesh& mesh){
    
    // itype v_start;
    // itype v_end;

    // n.get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range..

    // if(v_start == v_end) //no internal vertices..
    //     return;
    // itype v_range=v_end-v_start;
     
      leaf_VT vts;//(v_end-v_start,VT());

    n.get_VT(vts,dom,mesh);
    //   for(unsigned i=0;i<v_range;i++)
    // {
    //    itype real_v_id=v_start+i;

    //     VT &vt=vts[i];}
    
    }
    
    
    
    void Gradient::VV_relation(Node_V& n, Mesh& mesh, Spatial_Subdivision& division){
    
                /// if there are no vertices in the leaf we have nothing to do..
      if (n.is_leaf())
    {
        this->VV_relation_leaf(n,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->VV_relation(*n.get_son(i),mesh,division);
            }
        }
    }
    
    }
    
    void Gradient::VV_relation(Node_T &n, Box &n_dom, int level, Mesh &mesh, Spatial_Subdivision &division){
        
            if (n.is_leaf())
    {
        this->VV_relation_leaf(n,n_dom,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(n_dom,level,i);
            int son_level = level +1;
            if(n.get_son(i)!=NULL)
            {
                this->VV_relation(*n.get_son(i),son_dom,son_level,mesh,division);
            }
        }
    }
    }
    
    void Gradient::VV_relation_leaf(Node_V& n, Mesh& mesh){
    
//    if(!n.indexes_vertices())
//         return;

//     itype v_start = n.get_v_start();
//     itype v_end = n.get_v_end();
//     itype v_range = v_end - v_start;

    leaf_VV vvs;//(v_range,VV());
    //dvect v_aux(v_range,0.0);
    n.get_VV(vvs,mesh);
    // for(unsigned i=0;i<v_range;i++)
    // {
      
    //     itype real_v_id=v_start+i;
    //     Vertex &v = mesh.get_vertex(real_v_id);
    //     VV &vv = vvs[i];   
    // }
    }
    
    void Gradient::VV_relation_leaf(Node_T& n, Box &dom, Mesh& mesh){
    
    //  itype v_start;
    // itype v_end;

    // n.get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range..
    
    // if(v_start == v_end) //no internal vertices..
    //     return;
   
 //   leaf_VT vts(v_end-v_start,VT());
    
    leaf_VV vvs;//(v_end-v_start,VV());
    //dvect v_aux(v_range,0.0);
    n.get_VV(vvs,dom,mesh);
    
     
    // for(unsigned i=0;i<vvs.size();i++)
    // {
    //     itype real_v_id=v_start+i;
    //     Vertex &v = mesh.get_vertex(real_v_id);
    //     VV &vv = vvs[i];
    
    }
    
    }