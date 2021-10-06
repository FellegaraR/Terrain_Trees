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

#ifndef FORMANGRADIENTVECTOR_STATS_H
#define FORMANGRADIENTVECTOR_STATS_H

#include "forman_gradient_aux_structure.h"
#include "forman_gradient_caching_structure.h"
using namespace forman_aux_structures;

class FGV_cache_dangling_paths_stats
{
public:
    FGV_cache_dangling_paths_stats()
    {
        max_cache_items = 0;
        max_cache_leaves = 0;
        max_cache_internal_map_entry = 0;

        max_cache2_items = 0;
        max_cache2_leaves = 0;
        max_cache2_internal_map_entry = 0;

        max_dangling_desc2cells_path_item = 0;
        max_dangling_desc1cells_path_item = 0;
        max_dangling_asc2cells_path_item = 0;
        max_dangling_asc1cells_path_item = 0;

        max_dangling_desc2cells_path_leaves = 0;
        max_dangling_desc1cells_path_leaves = 0;
        max_dangling_asc2cells_path_leaves = 0;
        max_dangling_asc1cells_path_leaves = 0;

        max_queue_size = 0;

        max_local_MIG_arcs = max_local_MIG_minima = max_local_MIG_saddle = max_local_MIG_maxima = max_local_MIG_critica = 0;
    }

    inline void reset_stats_counter()
    {
        max_cache_items = 0;
        max_cache_leaves = 0;
        max_cache_internal_map_entry = 0;

        max_cache2_items = 0;
        max_cache2_leaves = 0;
        max_cache2_internal_map_entry = 0;

        max_dangling_desc2cells_path_item = 0;
        max_dangling_desc1cells_path_item = 0;
        max_dangling_asc2cells_path_item = 0;
        max_dangling_asc1cells_path_item = 0;

        max_dangling_desc2cells_path_leaves = 0;
        max_dangling_desc1cells_path_leaves = 0;
        max_dangling_asc2cells_path_leaves = 0;
        max_dangling_asc1cells_path_leaves = 0;

        max_queue_size = 0;

        max_local_MIG_arcs = max_local_MIG_minima = max_local_MIG_saddle = max_local_MIG_maxima = max_local_MIG_critica = 0;
    }

    inline void set_cache_stats(asc3_cache &cache)
    {
        set_cache_stats(cache.get_vtstar_cache());

        max_cache2_leaves = max(cache.get_vv_cache_size(),max_cache2_leaves);

        utype actual_entries_counter = 0;
        utype actual_items_counter = 0;

        for(vv_cache::mapIt it=cache.begin_vv_cache(); it!=cache.end_vv_cache(); ++it)
        {
            leaf_VV &vv = it->second;

            actual_entries_counter += vv.size();

            for(leaf_VV::iterator it=vv.begin(); it!=vv.end(); ++it)
                actual_items_counter += it->size();
        }

        max_cache2_internal_map_entry = max(max_cache2_internal_map_entry,actual_entries_counter);
        max_cache2_items = max(max_cache2_items,actual_items_counter);
    }

    inline void set_cache_stats(et_cache &cache)
    {
        max_cache_leaves = max(cache.get_actual_size(),max_cache_leaves);

        utype actual_entries_counter = 0;

        for(et_cache::mapIt it=cache.begin(); it!=cache.end(); ++it)
        {
            actual_entries_counter += it->second.size();
        }

        max_cache_internal_map_entry = max(max_cache_internal_map_entry,actual_entries_counter);
        max_cache_items = max(max_cache_items,actual_entries_counter);
    }

    inline void set_cache_stats(vtstar_cache &cache)
    {
        max_cache_leaves = max(cache.get_actual_size(),max_cache_leaves);

        utype actual_entries_counter = 0;

        for(vtstar_cache::mapIt it=cache.begin(); it!=cache.end(); ++it)
        {
            actual_entries_counter += it->second.size();
        }

        max_cache_internal_map_entry = max(max_cache_internal_map_entry,actual_entries_counter);
        max_cache_items = max(max_cache_items,actual_entries_counter);
    }

    inline void set_cache_stats(mig_cache &cache)
    {
        utype actual_entries_counter = 0;

        //VTSTAR CACHE
        set_cache_stats(cache.get_vtstar_cache());

        //ET CACHE
        max_cache2_leaves = max(cache.get_et_cache().get_actual_size(),max_cache2_leaves);
        for(et_cache::mapIt it=cache.begin_et_cache(); it!=cache.end_et_cache(); ++it)
        {
            leaf_ET &ef = it->second;

            actual_entries_counter += ef.size();
        }
        max_cache2_internal_map_entry = max(max_cache2_internal_map_entry,actual_entries_counter);
    }

    inline void set_2desc_dangling_paths_stats(leaves_2_desc_map &dangling)
    {
        max_dangling_desc2cells_path_leaves = max(max_dangling_desc2cells_path_leaves,static_cast<utype>(dangling.size()));

        utype count=0;

        for(leaves_2_desc_map::iterator it = dangling.begin(); it!=dangling.end(); ++it)
        {
            count += it->second.size();
        }

        max_dangling_desc2cells_path_item = max(max_dangling_desc2cells_path_item,count);
    }

    inline void set_1desc_dangling_paths_stats(leaves_1_desc_map &dangling)
    {
        max_dangling_desc1cells_path_leaves = max(max_dangling_desc1cells_path_leaves,static_cast<utype>(dangling.size()));

        utype count=0;

        for(leaves_1_desc_map::iterator it = dangling.begin(); it!=dangling.end(); ++it)
        {
            count += it->second.size();
        }

        max_dangling_desc1cells_path_item = max(max_dangling_desc1cells_path_item,count);
    }
    inline void set_2asc_dangling_paths_stats(leaves_2_asc_map &dangling)
    {
        max_dangling_asc2cells_path_leaves = max(max_dangling_asc2cells_path_leaves,static_cast<utype>(dangling.size()));

        utype count=0;

        for(leaves_2_asc_map::iterator it = dangling.begin(); it!=dangling.end(); ++it)
        {
            count += it->second.size();
        }

        max_dangling_asc2cells_path_item = max(max_dangling_asc2cells_path_item,count);
    }
    inline void set_1asc_dangling_paths_stats(leaves_asc_map &dangling)
    {
        max_dangling_asc1cells_path_leaves = max(max_dangling_asc1cells_path_leaves,static_cast<utype>(dangling.size()));

        utype count=0;

        for(leaves_asc_map::iterator it = dangling.begin(); it!=dangling.end(); ++it)
        {
            count += it->second.size();
        }

        max_dangling_asc1cells_path_item = max(max_dangling_asc1cells_path_item,count);
    }

    inline void set_all_dangling_paths(ig_paths &paths)
    {
        /// ascending 1 manifolds dangling paths
        max_dangling_asc1cells_path_leaves = max(max_dangling_asc1cells_path_leaves,paths.size_asc1cells());
        utype count=0;
        for(leaves_1_asc_mig_map::iterator it = paths.begin_asc1cells(); it!=paths.end_asc1cells(); ++it)
        {
            count += it->second.size();
        }
        max_dangling_asc1cells_path_item = max(max_dangling_asc1cells_path_item,count);

        /// descending 1 manifolds dangling paths
        max_dangling_desc1cells_path_leaves = max(max_dangling_desc1cells_path_leaves,paths.size_desc1cell());
        count=0;
        for(leaves_1_desc_mig_map::iterator it = paths.begin_desc1cell(); it!=paths.end_desc1cell(); ++it)
        {
            count += it->second.size();
        }
        max_dangling_desc1cells_path_item = max(max_dangling_desc1cells_path_item,count);
    }

    inline void print_stats()
    {
        if(max_cache_items > 0)
            cerr<<"-- -- --"<<endl;
        if(max_cache_items > 0)
        {
            cerr<<"cache: "<<max_cache_leaves<<" ";
            if(max_cache_internal_map_entry > 0)
                cerr<<max_cache_internal_map_entry<<" ";
            if(max_cache_items > 0)
                cerr<<max_cache_items;
            cerr<<endl;
        }

        if(max_cache2_items > 0)
        {
            cerr<<"cache2: "<<max_cache2_leaves<<" ";
            if(max_cache2_internal_map_entry > 0)
                cerr<<max_cache2_internal_map_entry<<" ";
            if(max_cache2_items > 0)
                cerr<<max_cache2_items;
            cerr<<endl;
        }

        if(max_dangling_desc2cells_path_item > 0)
        {
            cerr<<"desc2cells paths: "<<max_dangling_desc2cells_path_leaves<<" "<<max_dangling_desc2cells_path_item<<endl;
        }
        if(max_dangling_desc1cells_path_leaves > 0)
        {
            cerr<<"desc1cells paths: "<<max_dangling_desc1cells_path_leaves<<" "<<max_dangling_desc1cells_path_item<<endl;
        }

        if(max_dangling_asc2cells_path_item > 0)
        {
            cerr<<"asc2cells paths: "<<max_dangling_asc2cells_path_leaves<<" "<<max_dangling_asc2cells_path_item<<endl;
        }
        if(max_dangling_asc1cells_path_item > 0)
        {
            cerr<<"asc1cells paths: "<<max_dangling_asc1cells_path_leaves<<" "<<max_dangling_asc1cells_path_item<<endl;
        }

        if(max_local_MIG_critica > 0)
        {
            cerr<<"-- local MIG --"<<endl;
            cerr<<"maxima local critica number: "<<max_local_MIG_critica<<endl;
            cerr<<"max local min-saddle-max: "<<max_local_MIG_minima<<" "<<max_local_MIG_saddle<<" "<<max_local_MIG_maxima<<endl;
            cerr<<"max local arcs number: "<<max_local_MIG_arcs<<endl;
        }

        if(max_queue_size > 0)
            cerr<<"max queue size: "<<max_queue_size<<endl;

        if(max_cache_items > 0)
            cerr<<"-- -- --"<<endl;
    }

    inline void set_local_MIG_stats(IG &ig)
    {
        if(ig.get_arcs_number() > max_local_MIG_arcs)
            max_local_MIG_arcs = ig.get_arcs_number();
        if(ig.minima_number() > max_local_MIG_minima)
            max_local_MIG_minima = ig.minima_number();
        if(ig.saddle_number() > max_local_MIG_saddle)
            max_local_MIG_saddle = ig.saddle_number();
        if(ig.maxima_number() > max_local_MIG_maxima)
            max_local_MIG_maxima = ig.maxima_number();
        if(ig.critica_number() > max_local_MIG_critica)
            max_local_MIG_critica = ig.critica_number();
    }

private:
    utype max_cache_items; //all the item (max number of int in cache)
    utype max_cache_internal_map_entry;
    utype max_cache_leaves; //max number of leaves cached

    /// --- this 3 values below are only needed by the cache used for ascending 2 manifolds and mig--- ///
    utype max_cache2_items; //all the item (max number of int in cache)
    utype max_cache2_internal_map_entry;
    utype max_cache2_leaves; //max number of leaves cached

    utype max_dangling_desc2cells_path_item;
    utype max_dangling_desc2cells_path_leaves;

    utype max_dangling_desc1cells_path_item;
    utype max_dangling_desc1cells_path_leaves;

    utype max_dangling_asc2cells_path_item;
    utype max_dangling_asc2cells_path_leaves;

    utype max_dangling_asc1cells_path_item;
    utype max_dangling_asc1cells_path_leaves;

    /// needed to keep simplification stats
    utype max_local_MIG_minima, max_local_MIG_saddle, max_local_MIG_maxima;
    utype max_local_MIG_arcs, max_local_MIG_critica;

public:
    utype max_queue_size;
};

#endif // FORMANGRADIENTVECTOR_STATS_H
