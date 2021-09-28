/*
    This file is part of the Terrain Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)
               Federico Iuricich (federico.iuricich@gmail.com)

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

#include "forman_gradient_features_extractor.h"
#include "io/writer_morse.h"
#include "utilities/container_utilities.h"

void Forman_Gradient_Features_Extractor::extract_critical_clusters(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, map<short,ivect_set> &critical_simplices,
                                                                   Spatial_Subdivision &division, OpType operation, int cache_size, string mesh_name)
{
    vtstar_cache cache = vtstar_cache(cache_size);

    Timer time;

    critical_clusters cc = critical_clusters(mesh.get_vertices_num(),mesh.get_triangles_num(),
                                             critical_simplices[0].size(),critical_simplices[2].size());

    time.start();
    /// iterate over the minima..
    /// initialization phase
    for(ivect cV : critical_simplices[0])
    {
        cc.set_v_label(cV.front(),cc.get_counter_value());
        cc.add_vPerL(cc.get_counter_value(),cV.front());
        cc.increment_counter();
    }
    cc.reset_counter();
    ivect_set maxima = critical_simplices[2];
    this->init_critical_t_labels(n,mesh,gradient,maxima,division,cc);

    ivect_set saddles = critical_simplices[1];
    this->critical_clusters_extraction(n,mesh,gradient,saddles,division,n,cache,cc);

    for(int i=1; i<=mesh.get_vertices_num(); i++)
        cc.set_v_label(i,cc.get_vPerL_size(cc.get_v_label(i)));
    for(int i=1; i<=mesh.get_triangles_num(); i++)
        cc.set_t_label(i,cc.get_tPerL_size(cc.get_t_label(i)));

    time.stop();
    time.print_elapsed_time("[TIME] Extracting Critical Clusters: ");

    cerr<<"[STATS] cluster-min: "<<critical_simplices[0].size()-cc.get_cmin_counter_value()<<" ";
    cerr<<"cluster-max: "<<critical_simplices[2].size()-cc.get_cmax_counter_value()<<endl;

    if(operation == OUTPUT)
        Writer_Morse::write_critical_clusters(mesh_name,cc,mesh);
}

void Forman_Gradient_Features_Extractor::extract_critical_clusters(Node_T &n, Mesh &mesh, Forman_Gradient &gradient, map<short, ivect_set> &critical_simplices,
                                                                   Spatial_Subdivision &division, OpType operation, int cache_size, string mesh_name)
{
    Timer time;
    vtstar_cache cache = vtstar_cache(cache_size);
    critical_clusters cc = critical_clusters(mesh.get_vertices_num(),mesh.get_triangles_num(),
                                             critical_simplices[0].size(),critical_simplices[2].size());
    time.start();

    /// iterate over the minima..
    /// initialization phase
    for(ivect cV : critical_simplices[0])
    {
        cc.set_v_label(cV.front(),cc.get_counter_value());
        cc.add_vPerL(cc.get_counter_value(),cV.front());
        cc.increment_counter();
    }
    cc.reset_counter();
    ivect_set maxima = critical_simplices[2];
    this->init_critical_t_labels(n,mesh.get_domain(),0,mesh,gradient,maxima,division,cc);

    ivect_set saddles = critical_simplices[1];
    this->critical_clusters_extraction(n,mesh.get_domain(),0,mesh,gradient,saddles,division,n,cache,cc);

    for(int i=1; i<=mesh.get_vertices_num(); i++)
        cc.set_v_label(i,cc.get_vPerL_size(cc.get_v_label(i)));
    for(int i=1; i<=mesh.get_triangles_num(); i++)
        cc.set_t_label(i,cc.get_tPerL_size(cc.get_t_label(i)));

    time.stop();
    time.print_elapsed_time("[TIME] Extracting Critical Clusters: ");

    cerr<<"[STATS] cluster-min: "<<critical_simplices[0].size()-cc.get_cmin_counter_value()<<" ";
    cerr<<"cluster-max: "<<critical_simplices[2].size()-cc.get_cmax_counter_value()<<endl;

    if(operation == OUTPUT)
        Writer_Morse::write_critical_clusters(mesh_name,cc,mesh);
}

void Forman_Gradient_Features_Extractor::critical_clusters_extraction(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, ivect_set &saddles,
                                                                      Spatial_Subdivision &division, Node_V &root, vtstar_cache &cache, critical_clusters &cc)
{
    if (n.is_leaf())
    {
        critical_clusters_extraction_leaf(n,mesh,gradient,saddles,division,root,cache,cc);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->critical_clusters_extraction(*n.get_son(i),mesh,gradient,saddles,division,root,cache,cc);
            }
        }
    }
}

void Forman_Gradient_Features_Extractor::critical_clusters_extraction(Node_T &n, Box &n_dom, int n_level, Mesh &mesh, Forman_Gradient &gradient, ivect_set &saddles,
                                                                      Spatial_Subdivision &division, Node_T &root, vtstar_cache &cache, critical_clusters &cc)
{
    if (n.is_leaf())
    {
        critical_clusters_extraction_leaf(n,n_dom,mesh,gradient,saddles,division,root,cache,cc);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(n_dom,n_level,i);
            int son_level = n_level +1;
            if(n.get_son(i)!=NULL)
            {
                this->critical_clusters_extraction(*n.get_son(i),son_dom,son_level,mesh,gradient,saddles,division,root,cache,cc);
            }
        }
    }
}

void Forman_Gradient_Features_Extractor::critical_clusters_extraction_leaf(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, ivect_set &saddles,
                                                                           Spatial_Subdivision &division, Node_V &root, vtstar_cache &cache, critical_clusters &cc)
{
    if(!n.indexes_vertices())
        return;

    local_VTstar_ET local_rels;
    Forman_Gradient_Topological_Relations::get_VTstar_ET(local_rels,n,mesh,gradient);

    for(ivect_set::iterator c = saddles.begin(); c != saddles.end(); )
    {
        ivect cS = *c;

        if(!n.visited_vertex(cS.front())) // when I found the first critical simplex not indexed by the current node, I am done
            break;

        if(n.indexes_vertex(cS.back()))
        {
            clusterize_minima(*c,local_rels.get_VTstars(),cache,n,mesh,gradient,division,root,cc);
            clusterize_maxima(*c,local_rels.get_ETs(),gradient,cc);
            saddles.erase(c++); /// this invalidate the global structure encoding the critical simplices.. but makes efficient the computation
        }
        else
            ++c;
    }

    cache.insert(n.get_v_start()+n.get_v_end(),local_rels.get_VTstars());
}

void Forman_Gradient_Features_Extractor::critical_clusters_extraction_leaf(Node_T &n, Box &n_dom, Mesh &mesh, Forman_Gradient &gradient, ivect_set &saddles,
                                                                           Spatial_Subdivision &division, Node_T &root, vtstar_cache &cache, critical_clusters &cc)
{
    itype v_start;
    itype v_end;

    n.get_v_range(v_start,v_end,n_dom,mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;

    local_VTstar_ET local_rels;
    Forman_Gradient_Topological_Relations::get_VTstar_ET(local_rels,n,v_start,v_end,mesh,gradient);

    for(ivect_set::iterator c = saddles.begin(); c != saddles.end(); )
    {
        ivect cS = *c;

        if(cS.front() >= v_end) // when I found the first critical simplex not indexed by the current node, I am done
            break;

        if(n.indexes_vertex(v_start,v_end,cS.back()))
        {
            clusterize_minima(*c,local_rels.get_VTstars(),cache,n,v_start,v_end,mesh,gradient,division,root,cc);
            clusterize_maxima(*c,local_rels.get_ETs(),gradient,cc);
            saddles.erase(c++); /// this invalidate the global structure encoding the critical simplices.. but makes efficient the computation
        }
        else
            ++c;

    }

    cache.insert(v_start+v_end,local_rels.get_VTstars());
}

void Forman_Gradient_Features_Extractor::init_critical_t_labels(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, ivect_set &maxima,
                                                                Spatial_Subdivision &division, critical_clusters &cc)
{
    if (n.is_leaf())
    {
        if(!n.indexes_vertices())
            return;

        leaf_VT vts;
        n.get_VT(vts,mesh);

        for(ivect_set::iterator c = maxima.begin(); c != maxima.end(); )
        {
            ivect cS = *c;

            if(n.indexes_vertex(cS.back()))
            {
                itype t_id = Forman_Gradient_Topological_Relations::get_triangle_id(cS,vts[cS.back()-n.get_v_start()],mesh);
                cc.set_t_label(t_id,cc.get_counter_value());
                cc.add_tPerL(cc.get_counter_value(),t_id);
                cc.increment_counter();

                maxima.erase(c++); /// this invalidate the global structure encoding the critical simplices.. but makes efficient the computation
            }
            else
                ++c;
            if(!n.visited_vertex(cS.front())) // when I found the first critical simplex not indexed by the current node, I am done
                break;
        }
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->init_critical_t_labels(*n.get_son(i),mesh,gradient,maxima,division,cc);
            }
        }
    }
}

void Forman_Gradient_Features_Extractor::init_critical_t_labels(Node_T &n, Box &n_dom, int level, Mesh &mesh, Forman_Gradient &gradient, ivect_set &maxima,
                                                                Spatial_Subdivision &division, critical_clusters &cc)
{
    if (n.is_leaf())
    {
        itype v_start;
        itype v_end;

        n.get_v_range(v_start,v_end,n_dom,mesh); // we need to gather the vertices range..

        if(v_start == v_end) //no internal vertices..
            return;

        leaf_VT vts;
        n.get_VT(vts,v_start,v_end,mesh);

        for(ivect_set::iterator c = maxima.begin(); c != maxima.end(); )
        {
            ivect cS = *c;

            if(n.indexes_vertex(v_start,v_end,cS.back()))
            {
                itype t_id = Forman_Gradient_Topological_Relations::get_triangle_id(cS,vts[cS.back()-v_start],mesh);
                cc.set_t_label(t_id,cc.get_counter_value());
                cc.add_tPerL(cc.get_counter_value(),t_id);
                cc.increment_counter();

                maxima.erase(c++); /// this invalidate the global structure encoding the critical simplices.. but makes efficient the computation
            }
            else
                ++c;
            if(cS.front() >= v_end) // when I found the first critical simplex not indexed by the current node, I am done
                break;
        }
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(n_dom,level,i);
            int son_level = level +1;
            if(n.get_son(i)!=NULL)
            {
                this->init_critical_t_labels(*n.get_son(i),son_dom,son_level,mesh,gradient,maxima,division,cc);
            }
        }
    }
}

void Forman_Gradient_Features_Extractor::clusterize_minima(const ivect &saddle, leaf_VTstar &vtstars, vtstar_cache &cache, Node_V &n, Mesh &mesh,
                                                           Forman_Gradient &gradient, Spatial_Subdivision &division,
                                                           Node_V &root, critical_clusters &cc)
{
    itype t1_id = Forman_Gradient_Topological_Relations::get_VTstar(n,saddle.front(),vtstars,cache,root,division,mesh,gradient);
    itype t2_id = Forman_Gradient_Topological_Relations::get_VTstar(n,saddle.back(),vtstars,cache,root,division,mesh,gradient);

    if(gradient.is_vertex_critical(saddle.front(),t1_id,mesh) && gradient.is_vertex_critical(saddle.back(),t2_id,mesh))
    {
        if(cc.get_v_label(saddle.back()) != cc.get_v_label(saddle.front()))
        {

            itype first, second;
            if(cc.get_vPerL_size(cc.get_v_label(saddle.front())) > cc.get_vPerL_size(cc.get_v_label((saddle.back()))))
            {
                first = saddle.front();
                second = saddle.back();
            }
            else{
                second = saddle.front();
                first = saddle.back();
            }

            int label=cc.get_v_label(second);
            for(auto v : cc.get_vPerL(label)){
                cc.set_v_label(v,cc.get_v_label(first));
            }

            cc.add_vPerL(cc.get_v_label(first),cc.get_vPerL(label));
            cc.clear_vPerL(label);
            cc.increment_cmin_counter();
        }
    }
}

void Forman_Gradient_Features_Extractor::clusterize_minima(const ivect &saddle, leaf_VTstar &vtstars, vtstar_cache &cache, Node_T &n, itype v_start, itype v_end, Mesh &mesh, Forman_Gradient &gradient,
                                                           Spatial_Subdivision &division, Node_T &root, critical_clusters &cc)
{
    itype t1_id = Forman_Gradient_Topological_Relations::get_VTstar(n,v_start,v_end,saddle.front(),vtstars,cache,root,division,mesh,gradient);
    itype t2_id = Forman_Gradient_Topological_Relations::get_VTstar(n,v_start,v_end,saddle.back(),vtstars,cache,root,division,mesh,gradient);

    if(gradient.is_vertex_critical(saddle.front(),t1_id,mesh) && gradient.is_vertex_critical(saddle.back(),t2_id,mesh))
    {
        if(cc.get_v_label(saddle.back()) != cc.get_v_label(saddle.front()))
        {

            itype first, second;
            if(cc.get_vPerL_size(cc.get_v_label(saddle.front())) > cc.get_vPerL_size(cc.get_v_label((saddle.back()))))
            {
                first = saddle.front();
                second = saddle.back();
            }
            else{
                second = saddle.front();
                first = saddle.back();
            }

            int label=cc.get_v_label(second);
            for(auto v : cc.get_vPerL(label)){
                cc.set_v_label(v,cc.get_v_label(first));
            }

            cc.add_vPerL(cc.get_v_label(first),cc.get_vPerL(label));
            cc.clear_vPerL(label);
            cc.increment_cmin_counter();
        }
    }
}

void Forman_Gradient_Features_Extractor::clusterize_maxima(const ivect &saddle, leaf_ET &ets, Forman_Gradient &gradient, critical_clusters &cc)
{
    ET &et = ets[saddle];

    if(et.second!=-1 && gradient.is_triangle_critical(et.first) && gradient.is_triangle_critical(et.second))
    {
//        // for debug only //
//        this->saddles_around_maxima.insert(make_pair(saddle,et));
        //
        if(cc.get_t_label(et.first) != cc.get_t_label(et.second))
        {
            itype first, second;

            if(cc.get_tPerL_size(cc.get_t_label(et.first)) > cc.get_tPerL_size(cc.get_t_label(et.second)))
            {
                first = et.first;
                second = et.second;
            }
            else
            {
                second = et.first;
                first = et.second;
            }

            int label=cc.get_t_label(second);
            for(auto t : cc.get_tPerL(label)){
                cc.set_t_label(t,cc.get_t_label(first));
            }
            cc.add_tPerL(cc.get_t_label(first),cc.get_tPerL(label));
            cc.clear_tPerL(label);

            cc.increment_cmax_counter();
        }
    }
}
