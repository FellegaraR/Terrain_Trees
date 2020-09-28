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

void Forman_Gradient_Features_Extractor::extract_descending_2cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division, Node_V &root,
                                                                   OpType operation, int cache_size)
{
    /// reset the statistical counters
    reset_extraction_critica_counters();

    leaves_2_desc_map dangling_paths = leaves_2_desc_map();
    et_cache cache = et_cache(cache_size);
    this->descending_2cells_extraction(n,mesh,gradient,division,root,operation,dangling_paths,cache);

    if(dangling_paths.size() != 0)
    {
        cout<<"[ERROR] the procedure to not visit fully the dangling paths"<<endl;
        cout<<dangling_paths.size()<<endl;
        for(leaves_2_desc_map::iterator it=dangling_paths.begin(); it!=dangling_paths.end(); ++it)
            cout<<it->first<<" -- "<<it->second.size()<<endl;
        int a; cin>>a;
    }

    cache.reset();
    dangling_paths.clear();
    cout<<"maxima visited: "<<found_max<<endl;
    found_max = 0;
}

void Forman_Gradient_Features_Extractor::descending_2cells_extraction(Node_V &n, Mesh &mesh,  Forman_Gradient &gradient, Spatial_Subdivision &division, Node_V &root,
                                                                      OpType operation, leaves_2_desc_map &dangling_paths, et_cache &cache)
{
    if (n.is_leaf())
    {
        this->descending_2cells_extraction_leaf(n,mesh,gradient,root,division,operation,dangling_paths,cache);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->descending_2cells_extraction(*n.get_son(i),mesh,gradient,division,root,operation,dangling_paths,cache);
            }
        }
    }
}

void Forman_Gradient_Features_Extractor::descending_2cells_extraction_leaf(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Node_V &root, Spatial_Subdivision &division,
                                                                           OpType operation, leaves_2_desc_map &dangling_paths, et_cache &cache)
{
    /// if there are no vertices in the leaf we have nothing to do..
    if(!n.indexes_vertices())
        return;

    Timer time;
    if(operation == TIME_VERBOSE)
        time.start();
    leaf_ET local_ef;
    n.get_ET(local_ef,mesh);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        rels_gathering += time.get_elapsed_time();
        time.start();
    }

    get_new_descending_2cells(n,mesh,gradient,local_ef,root,division,operation,dangling_paths,cache);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        new_paths += time.get_elapsed_time();
        time.start();
    }

    get_dangling_descending_2cells(n,mesh,gradient,local_ef,root,division,operation,dangling_paths,cache);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        u_paths += time.get_elapsed_time();
        time.start();
    }

    //add the current vt to cache (in any leaves except the latest)
    if(mesh.get_vertices_num() > n.get_v_end())
        cache.insert(n.get_v_end()+n.get_v_start(),local_ef);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        cache_t += time.get_elapsed_time();
    }

    //get stats
    if(operation == OUTPUT)
    {
        stats.set_cache_stats(cache);
        stats.set_2desc_dangling_paths_stats(dangling_paths);
    }
}

void Forman_Gradient_Features_Extractor::get_new_descending_2cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, leaf_ET &local_ef, Node_V &root, Spatial_Subdivision &division,
                                                                   OpType operation, leaves_2_desc_map &dangling_paths, et_cache &cache)
{
    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        // we start a visit only when we visit the leaf block indexing the maximum vertex of the maximum
        if(n.indexes_vertex(mesh.get_triangle(*t_id).maxindex()) && gradient.is_triangle_critical(*t_id))
        {
            found_max++;
            get_one_descending_2cells(n, *t_id, mesh, gradient, *t_id, local_ef, root, division,operation,dangling_paths,cache);
        }
    }
}

void Forman_Gradient_Features_Extractor::get_dangling_descending_2cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, leaf_ET &local_et, Node_V &root, Spatial_Subdivision &division,
                                                                        OpType operation, leaves_2_desc_map &dangling_paths, et_cache &cache)
{

    //quindi cerco path iniziati in altre foglie che continuano in quella corrente
    int key = n.get_v_start()+n.get_v_end();
    leaves_2_desc_map::iterator it = dangling_paths.find(key);
    if(it != dangling_paths.end())
    {
        map_2paths &tri = it->second;
        map_2paths::iterator it_t;

        for(it_t = tri.begin(); it_t != tri.end();++it_t)
        {
            get_one_descending_2cells(n, it_t->first, mesh, gradient, it_t->second, local_et, root,division,operation,dangling_paths,cache);
        }

        dangling_paths.erase(it);
    }
}

void Forman_Gradient_Features_Extractor::get_one_descending_2cells(Node_V &n, itype t_id, Mesh &mesh, Forman_Gradient &gradient, itype label, leaf_ET &local_ef,
                                                                   Node_V &root, Spatial_Subdivision &division, OpType operation,
                                                                   leaves_2_desc_map &dangling_paths, et_cache &cache)
{
    if(operation == OUTPUT)
        segmentation[t_id-1] = label;

    iqueue coda;
    coda.push(t_id);
    itype current_tri;

    while(!coda.empty())
    {
        current_tri = coda.front();
        coda.pop();

        Triangle& t = mesh.get_triangle(current_tri);

        //if I have all the vertices visited, then I have all the topological relation locally defined in n or (hopefully) in cache
        if(n.has_all_vertices_visited(t))
        {
            for(int f=0; f<t.vertices_num(); f++)
            {
                ivect edge;
                t.TE(f,edge);

                pair<itype,itype> et = Forman_Gradient_Topological_Relations::get_ET(n,edge,local_ef,cache,root,division,mesh);

                /// with the new encoding I don't need to check that the face is paired in the adjacent triangle
                /// I already know it.. (I should)
                itype next = (et.first == current_tri) ? et.second : et.first;

                if(next != -1)
                {
                    TriGradient t_grad = gradient.convert_compressed_to_expand(next);
                    Triangle& t_next = mesh.get_triangle(next);

                    if(t_grad.get_edge_pair(t_next.vertex_index(edge[0]),t_next.vertex_index(edge[1])) == 3)
                    {
                        coda.push(next);
                        if(operation == OUTPUT)
                        {
                            segmentation[next-1] = label; //come etichetta metto l'id del triangolo critico da cui sono partito
                        }
                    }

                }
            }
        }
        else
        {            
            Dangling_Paths_Handler::add_dangling_path(root,division,dangling_paths,current_tri,t,label);
        }

        if(operation == OUTPUT && stats.max_queue_size < coda.size())
            stats.max_queue_size = coda.size();
    }
}
