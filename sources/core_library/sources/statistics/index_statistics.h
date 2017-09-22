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

#ifndef _INDEXSTATISTICS_H
#define	_INDEXSTATISTICS_H

#include <vector>
#include "basic_types/basic_wrappers.h"

using namespace std;

///A class representing a container used to store the statistics obtained from a spatial index
class IndexStatistics{
public:
    ///A constructor method
    IndexStatistics()
    {
        numNode=0;
        numFullLeaf=numEmptyLeaf=0;

        minTreeDepth=-1;
        avgTreeDepth=maxTreeDepth=0;
        min_vertices_per_leaf=-1;
        avg_vertices_per_leaf=max_vertices_per_leaf=0;

        min_partially_indexed_tri=-1;
        avg_partially_indexed_tri=max_partially_indexed_tri=0;

        min_overlapping_tri=-1;
        avg_overlapping_tri=max_overlapping_tri=0;

        min_completely_indexed_tri=-1;
        avg_completely_indexed_tri=max_completely_indexed_tri=0;

        numTin1Leaf = numTin2Leaf = numTin3Leaf = numTin4Leaf = numTinMoreLeaf = 0;
        min_leaves_for_tri=-1;
        avg_leaves_for_tri=max_leaves_for_tri=0;
        avg_weighted_leaves_for_tri = 0;

        t_list_length = 0;
        real_t_list_length = 0;
    }

    ///A public variable representing the number of tree nodes
    int numNode;
    ///A public variable representing the number of leaf containing significant information
    int numFullLeaf;
    ///A public variable representing the number of empty leaf
    int numEmptyLeaf;
    ///A public variable representing the minimum tree depth
    int minTreeDepth;
    ///A public variable representing the maximum tree depth
    int maxTreeDepth;
    ///A public variable representing the average tree depth
    coord_type avgTreeDepth;
    ///A public variable representing the minimum number of vertex per full leaf node
    int min_vertices_per_leaf;
    ///A public variable representing the m aximum number of vertex per full leaf node
    int max_vertices_per_leaf;
    ///A public variable representing the average number of vertex per full leaf node
    coord_type avg_vertices_per_leaf;
    ///A public variable representing the minimum number of triangles partially indexed by a leaf node
    int min_partially_indexed_tri;
    ///A public variable representing the maximum number of triangles partially indexed by a leaf node
    int max_partially_indexed_tri;
    ///A public variable representing the average number of triangles partially indexed by a leaf node
    coord_type avg_partially_indexed_tri;
    ///A public variable representing the minimum number of triangles overlapping-only a leaf node
    int min_overlapping_tri;
    ///A public variable representing the maximum number of triangles overlapping-only a leaf node
    int max_overlapping_tri;
    ///A public variable representing the average number of triangles overlapping-only a leaf node
    coord_type avg_overlapping_tri;
    ///A public variable representing the minimum number of triangles completely indexed by a leaf node
    int min_completely_indexed_tri;
    ///A public variable representing the maximum number of triangles completely indexed by a leaf node
    int max_completely_indexed_tri;
    ///A public variable representing the average number of triangles completely indexed by a leaf node
    coord_type avg_completely_indexed_tri;
    ///A public array, with an entry for each triangle, containing the number of leaf indexing each triangle
    ivect num_leaves_for_tri;
    ///A public variable representing the number of triangles indexed in exactly one leaf
    int numTin1Leaf;
    ///A public variable representing the number of triangles indexed in exactly two leaves
    int numTin2Leaf;
    ///A public variable representing the number of triangles indexed in exactly three leaves
    int numTin3Leaf;
    ///A public variable representing the number of triangles indexed in exactly four leaves
    int numTin4Leaf;
    ///A public variable representing the number of triangles indexed in more than four leaves
    int numTinMoreLeaf;
    ///A public variable representing the minimum number of leaves indexing a single triangle
    int min_leaves_for_tri;
    ///A public variable representing the maximum number of leaves indexing a single triangle
    int max_leaves_for_tri;
    ///A public variable representing the average number of leaves indexing a single triangle (chi)
    coord_type avg_leaves_for_tri;
    ///A public variable representing the average weighted number of leaves indexing a single triangle (chi_w)
    coord_type avg_weighted_leaves_for_tri;
    ///A public variable representing the summation of the compressed triangles arrays
    int t_list_length;
    ///A public variable representing the summation of the un-compressed triangles arrays
    int real_t_list_length;
};

#endif	/* _INDEXSTATISTICS_H */
