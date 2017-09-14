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

#ifndef TOPOLOGICALQUERIES_H
#define TOPOLOGICALQUERIES_H

#include <set>
#include <map>
#include <boost/dynamic_bitset.hpp>

#include "basic_types/vertex.h"
#include "basic_types/triangle.h"
#include "basic_types/mesh.h"
#include "basic_types/box.h"
#include "terrain_trees/node_v.h"
#include "terrain_trees/node_t.h"
#include "terrain_trees/spatial_subdivision.h"
#include "geometry/geometry_wrapper.h"

using namespace std;

/**
 * @brief The Topological_Queries class provides an interface for executing topological queries on the Terrain Trees
 * NOTA: of this class are documented only the public procedures.
 */
class Topological_Queries
{
public:
    /**
     * @brief A constructor method
     */
    Topological_Queries() { }
    /**
     * @brief A copy-constructor method
     */
    Topological_Queries(const Topological_Queries&) {}
    /**
     * @brief A destructor method
     */
    virtual ~Topological_Queries() {}
    ///A public method that excutes windowed VT queries, reading the boxes from file
    /*!
     * This method prints the results on standard output
     *
     * \param n a N& argument, representing the actual node to visit
     * \param dom a Box& argument, representing the node domain
     * \param mesh a Mesh& argument, representing the current mesh
     * \param division a D& argument, representing the tree subdivision type
     * \param query_path a string argument, representing the file path of the query input
     * \param reindexed a boolean, true if the index and the mesh are spatially reordered
     */
    template<class N> void windowed_VT(N &n, Box &dom, Mesh &mesh, Spatial_Subdivision &division, string query_path, bool reindexed);
    ///A public method that excutes windowed TT queries, reading the boxes from file
    /*!
     * This method prints the results on standard output
     *
     * \param n a N& argument, representing the actual node to visit
     * \param dom a Box& argument, representing the node domain
     * \param mesh a Mesh& argument, representing the current mesh
     * \param division a D& argument, representing the tree subdivision type
     * \param query_path a string argument, representing the file path of the query input
     */
    template<class N> void windowed_TT(N &n, Box &dom, Mesh &mesh, Spatial_Subdivision &division, string query_path);
    ///A public method that excutes batched VT queries
    /*!
     * This method prints the results on standard output
     *
     * \param n a N& argument, representing the actual node to visit
     * \param dom a Box& argument, representing the node domain
     * \param mesh a Mesh& argument, representing the current mesh
     * \param division a D& argument, representing the tree subdivision type
     * \param reindexed a boolean, true if the index and the mesh are spatially reordered
     */
    template<class N> void batched_VT(N &n, Box &dom, Mesh &mesh, Spatial_Subdivision &division, bool reindex);
    ///A public method that excutes batched TT queries
    /*!
     * This method prints the results on standard output
     *
     * \param n a N& argument, representing the actual node to visit
     * \param dom a Box& argument, representing the node domain
     * \param mesh a Mesh& argument, representing the current mesh
     * \param division a D& argument, representing the tree subdivision type
     * \param reindexed a boolean, true if the index and the mesh are spatially reordered
     */
    template<class N> void batched_TT(N &n, Mesh &mesh, Spatial_Subdivision &division);

private:
    // windowed VT - auxiliary functions
    template<class D> void windowed_VT(Node_T &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, wVT &vt);
    template<class N, class D> void windowed_VT_no_reindex(N &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, wVT &vt);
    template<class D> void windowed_VT(Node_V &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, wVT &vt);
    void windowed_VT_Leaf(Node_T& n, Box &dom, Box &b, Mesh& mesh, wVT &vt);
    void windowed_VT_Leaf(Node_V& n, Box &b, Mesh& mesh, wVT &vt);
    template<class N> void windowed_VT_Leaf_no_reindex(N& n, Box &dom, Box &b, Mesh& mesh, wVT &vt);
    void update_resulting_VT(itype v, itype t, wVT &vt);

    // windowed TT - auxiliary functions
    template<class N, class D> void windowed_TT(N &n, Box &dom, int level, Box &b, Mesh &mesh, D &division, wTT &tt,
                                                boost::dynamic_bitset<> &checked_tri);
    template<class N> void windowed_TT_Leaf_test(N& n, Box &b, Mesh& mesh, wTT &tt, boost::dynamic_bitset<> &checked_tri);
    template<class N> void windowed_TT_Leaf_add(N& n, Mesh& mesh, wTT &tt, boost::dynamic_bitset<> &checked_tri);
    void add_edges(itype t_id, vector<edge_triangle_tuple> &tuples, Mesh &mesh, wTT::const_iterator &iter, wTT &tt);
    void pair_adjacent_triangles(vector<edge_triangle_tuple> &tuples, Mesh &mesh, wTT &tt);
    void update_resulting_TT(int pos, itype t1, itype t2, wTT &tt);
    void init_TT_entry(itype t1, utype size, wTT &tt);
    void finalize_TT_Leaf(vector<edge_triangle_tuple> &tuples, wTT &tt, Mesh &mesh);

    // batched VT - auxiliary functions
    template<class N, class D> void batched_VT_visit(N &n, Box &dom, int level, Mesh &mesh, D &division, bool stats, utype &max_entries);
    template<class N, class D> void batched_VT_no_reindex(N &n, Box &dom, int level, Mesh &mesh, D &division, bool stats, utype &max_entries);
    void batched_VT_leaf(Node_T &n, Box &dom, Mesh &mesh, bool stats, utype &max_entries);
    void batched_VT_leaf(Node_V &n, Box &, Mesh &mesh, bool stats, utype &max_entries);
    void batched_VT_no_reindex_leaf(Node_T &n, Box &dom, Mesh &mesh, bool stats, utype &max_entries);
    void batched_VT_no_reindex_leaf(Node_V &n, Box &dom, Mesh &mesh, bool stats, utype &max_entries);

    // batched TT - auxiliary functions
    template<class N, class D> void batched_TT_visit(N &n, Mesh &mesh, D &division, vector<ivect > &tt, bool stats, utype &max_entries);
    template<class N> void batched_TT_leaf(N &n, Mesh &mesh, vector<ivect > &tt, bool stats, utype &max_entries);
};

#include "topological_queries_windowed.h"
#include "topological_queries_batched.h"

#endif // TOPOLOGICALQUERIES_H
