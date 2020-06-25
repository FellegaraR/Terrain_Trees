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

#ifndef FORMAN_GV_CACHING_STRUCTURE_H
#define FORMAN_GV_CACHING_STRUCTURE_H

#include "forman_gradient.h"
#include "basic_types/lru_cache.h"

namespace forman_aux_structures
{
typedef LRU_Cache<utype,leaf_VT> vt_cache;
typedef LRU_Cache<utype,leaf_VV> vv_cache;
typedef LRU_Cache<utype,leaf_VTstar> vtstar_cache;
typedef LRU_Cache<utype,leaf_ETstar> etstar_cache;
typedef LRU_Cache<utype,leaf_ET> et_cache;

class asc3_cache
{
protected:
    vv_cache vv_lru_cache;
    vtstar_cache vtstar_lru_cache;

public:
    asc3_cache(utype size)
    {
        vv_lru_cache.setMaxSize(size);
        vtstar_lru_cache.setMaxSize(size);
    }

    inline vv_cache::mapIt add_vv_lists(utype key, leaf_VV &lists) { return vv_lru_cache.insert(key,lists); }
    inline vv_cache::mapIt find_vv_lists(utype key) { return vv_lru_cache.find(key); }
    inline vv_cache::mapIt end_vv_cache() { return vv_lru_cache.end(); }
    inline vv_cache::mapIt begin_vv_cache() { return vv_lru_cache.begin(); }

    inline vtstar_cache::mapIt add_vtstar_lists(utype key, leaf_VTstar &lists) { return vtstar_lru_cache.insert(key,lists); }
    inline vtstar_cache::mapIt find_vtstar_lists(utype key) { return vtstar_lru_cache.find(key); }
    inline vtstar_cache::mapIt end_vtstar_cache() { return vtstar_lru_cache.end(); }
    inline vtstar_cache::mapIt begin_vtstar_cache() { return vtstar_lru_cache.begin(); }

    vtstar_cache& get_vtstar_cache() { return vtstar_lru_cache; }
    utype get_vv_cache_size() { return vv_lru_cache.get_actual_size(); }

    inline void reset()
    {
        vv_lru_cache.reset();
        vtstar_lru_cache.reset();
    }
    inline void addItem(utype key, asc2rels &rels)
    {
        vv_lru_cache.insert(key,rels.vvs);
        vtstar_lru_cache.insert(key,rels.vtstars);
    }

};

class mig_cache
{
private:
    vtstar_cache vtstar_lru_cache;
    et_cache et_lru_cache;

public:
    mig_cache(utype size)
    {
        vtstar_lru_cache.setMaxSize(size);
        et_lru_cache.setMaxSize(size);
    }

    inline vtstar_cache::mapIt add_vtstar_lists(itype key, leaf_VTstar &lists) { return vtstar_lru_cache.insert(key,lists); }
    inline vtstar_cache::mapIt find_vtstar_lists(itype key) { return vtstar_lru_cache.find(key); }
    inline vtstar_cache::mapIt end_vtstar_cache() { return vtstar_lru_cache.end(); }
    inline vtstar_cache::mapIt begin_vtstar_cache() { return vtstar_lru_cache.begin(); }

    inline et_cache::mapIt add_ef_lists(itype key, leaf_ET &lists) { return et_lru_cache.insert(key,lists); }
    inline et_cache::mapIt find_ef_lists(itype key) { return et_lru_cache.find(key); }
    inline et_cache::mapIt end_et_cache() { return et_lru_cache.end(); }
    inline et_cache::mapIt begin_et_cache() { return et_lru_cache.begin(); }

    vtstar_cache& get_vtstar_cache() { return vtstar_lru_cache; }
    et_cache& get_et_cache() { return et_lru_cache; }

    inline void reset()
    {
        vtstar_lru_cache.reset();
        et_lru_cache.reset();
    }
    inline void addItem(itype key, local_VTstar_ET &rels)
    {
        vtstar_lru_cache.insert(key,rels.get_VTstars());
        et_lru_cache.insert(key,rels.get_ETs());
    }
};

}

#endif // FORMAN_GV_CACHING_STRUCTURE_H
