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

#include "statistics.h"

void Statistics::compute_leaf_statistics(Node_T &n, Box& dom, Mesh& mesh, bool)
{
    int num_t_completely = 0;
    int num_t_partially = 0;
    int num_t_overlapping = 0;

    this->indexStats.t_list_length += n.get_t_array_size();
    this->indexStats.real_t_list_length += n.get_real_t_array_size();

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& runIt = itPair.first;
        Triangle& tri = mesh.get_triangle(*runIt);
        if(n.completely_indexes_triangle_vertices_dom(tri,dom,mesh))
            num_t_completely++;
        else if(n.indexes_triangle_vertices_dom(tri,dom,mesh))
            num_t_partially++;
        else
            num_t_overlapping++;
        this->indexStats.num_leaves_for_tri[*runIt-1]++;
    }

    if((num_t_completely + num_t_overlapping + num_t_partially) > 0)
    {
        this->indexStats.numFullLeaf++;
        set_leaf_triangles_stats(num_t_completely,num_t_partially,num_t_overlapping);
    }
    else
    {
        this->indexStats.numEmptyLeaf++;
    }
}

void Statistics::compute_leaf_statistics(Node_V &n, Box& dom, Mesh& mesh, bool reindex)
{
    int num_t_completely = 0;
    int num_t_partially = 0;
    int num_t_overlapping = 0;

    this->indexStats.t_list_length += n.get_t_array_size();
    this->indexStats.real_t_list_length += n.get_real_t_array_size();

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& runIt = itPair.first;
        Triangle& tri = mesh.get_triangle(*runIt);
        if((reindex && n.completely_indexes_triangle_vertices(tri)) || n.completely_indexes_triangle_vertices_dom(tri,dom,mesh))
            num_t_completely++;
        else if((reindex && n.indexes_triangle_vertices(tri)) || n.indexes_triangle_vertices_dom(tri,dom,mesh))
            num_t_partially++;
        else
            num_t_overlapping++;
        this->indexStats.num_leaves_for_tri[*runIt-1]++;
    }


    if((num_t_completely + num_t_overlapping + num_t_partially) > 0)
    {
        this->indexStats.numFullLeaf++;
        set_leaf_vertices_stats(n.get_real_v_array_size());
        set_leaf_triangles_stats(num_t_completely,num_t_partially,num_t_overlapping);
    }
    else
    {        
        this->indexStats.numEmptyLeaf++;
    }
}

void Statistics::set_leaf_vertices_stats(int num_vertex)
{
    if(this->indexStats.min_vertices_per_leaf==-1 || this->indexStats.min_vertices_per_leaf > num_vertex)
        this->indexStats.min_vertices_per_leaf = num_vertex;
    if(this->indexStats.max_vertices_per_leaf < num_vertex)
        this->indexStats.max_vertices_per_leaf = num_vertex;
    this->indexStats.avg_vertices_per_leaf += num_vertex;
}

void Statistics::set_leaf_triangles_stats(int num_t_completely, int num_t_partially, int num_t_overlapping)
{
    if(this->indexStats.min_completely_indexed_tri==-1 || this->indexStats.min_completely_indexed_tri > num_t_completely)
        this->indexStats.min_completely_indexed_tri =  num_t_completely;
    if(this->indexStats.max_completely_indexed_tri < num_t_completely)
        this->indexStats.max_completely_indexed_tri = num_t_completely;
    this->indexStats.avg_completely_indexed_tri += num_t_completely;

    if(this->indexStats.min_partially_indexed_tri==-1 || this->indexStats.min_partially_indexed_tri > num_t_partially)
        this->indexStats.min_partially_indexed_tri = num_t_partially;
    if(this->indexStats.max_partially_indexed_tri < num_t_partially)
        this->indexStats.max_partially_indexed_tri = num_t_partially;
    this->indexStats.avg_partially_indexed_tri += num_t_partially;

    if(this->indexStats.min_overlapping_tri==-1 || this->indexStats.min_overlapping_tri > num_t_overlapping)
        this->indexStats.min_overlapping_tri = num_t_overlapping;
    if(this->indexStats.max_overlapping_tri < num_t_overlapping)
        this->indexStats.max_overlapping_tri = num_t_overlapping;
    this->indexStats.avg_overlapping_tri += num_t_overlapping;
}

void Statistics::calc_remaining_index_statistics()
{
    for(ivect_iter iter=this->indexStats.num_leaves_for_tri.begin(); iter!=this->indexStats.num_leaves_for_tri.end(); ++iter)
    {
        if(*iter == 1)
            this->indexStats.numTin1Leaf++;
        else if(*iter == 2)
            this->indexStats.numTin2Leaf++;
        else if(*iter == 3)
            this->indexStats.numTin3Leaf++;
        else if(*iter == 4)
            this->indexStats.numTin4Leaf++;
        else
            this->indexStats.numTinMoreLeaf++;

        if(this->indexStats.min_leaves_for_tri == -1 || this->indexStats.min_leaves_for_tri > *iter)
            this->indexStats.min_leaves_for_tri = *iter;
        if(this->indexStats.max_leaves_for_tri < *iter)
            this->indexStats.max_leaves_for_tri = *iter;

        this->indexStats.avg_leaves_for_tri += *iter;

        if(*iter!=1)
            this->indexStats.avg_weighted_leaves_for_tri += *iter;
    }

    if(this->indexStats.numNode > 0)
        this->indexStats.avgTreeDepth /= (this->indexStats.numEmptyLeaf + this->indexStats.numFullLeaf);
    if(this->indexStats.numFullLeaf > 0)
    {
        this->indexStats.avg_vertices_per_leaf /= this->indexStats.numFullLeaf;
        this->indexStats.avg_partially_indexed_tri /= this->indexStats.numFullLeaf;

        if(this->indexStats.avg_overlapping_tri != 0)
            this->indexStats.avg_overlapping_tri /= this->indexStats.numFullLeaf;
        if(this->indexStats.avg_completely_indexed_tri != 0)
            this->indexStats.avg_completely_indexed_tri /= this->indexStats.numFullLeaf;
    }
    if(this->indexStats.num_leaves_for_tri.size() > 0)
    {
        this->indexStats.avg_leaves_for_tri /= this->indexStats.num_leaves_for_tri.size();
        this->indexStats.avg_weighted_leaves_for_tri /= (this->indexStats.numTin2Leaf+this->indexStats.numTin3Leaf+this->indexStats.numTin4Leaf+this->indexStats.numTinMoreLeaf);
    }
}

void Statistics::check_inconsistencies()
{
    if(this->indexStats.numEmptyLeaf == 0)
        this->indexStats.min_overlapping_tri = 0;
    if(this->indexStats.numFullLeaf == 0)
    {
        this->indexStats.min_vertices_per_leaf = 0;
        this->indexStats.min_partially_indexed_tri = 0;
    }
    if(this->indexStats.numNode == 0)
        this->indexStats.minTreeDepth = 0;
    if(this->indexStats.num_leaves_for_tri.size() == 0)
        this->indexStats.avg_leaves_for_tri = 0;

    if(this->indexStats.min_vertices_per_leaf==-1)
        this->indexStats.min_vertices_per_leaf=0;
    if(this->indexStats.min_overlapping_tri==-1)
        this->indexStats.min_overlapping_tri=0;

    return;
}

int Statistics::compute_queries_statistics(QueryStatistics &qS)
{
    int hit_ratio=0;

    if(qS.triangles.size() > 0)
        hit_ratio++;

    //triangles statistics
    if(this->fullQueryStats.mintri > qS.triangles.size())
        this->fullQueryStats.mintri = qS.triangles.size();
    if(this->fullQueryStats.maxtri < qS.triangles.size())
        this->fullQueryStats.maxtri = qS.triangles.size();
    this->fullQueryStats.avgtri += qS.triangles.size();

    //nodes statistics
    if(this->fullQueryStats.minNode > qS.numNode)
        this->fullQueryStats.minNode = qS.numNode;
    if(this->fullQueryStats.maxNode < qS.numNode)
        this->fullQueryStats.maxNode = qS.numNode;
    this->fullQueryStats.avgNode += qS.numNode;

    //leaves statistics
    if(this->fullQueryStats.minLeaf > qS.numLeaf)
        this->fullQueryStats.minLeaf = qS.numLeaf;
    if(this->fullQueryStats.maxLeaf < qS.numLeaf)
        this->fullQueryStats.maxLeaf = qS.numLeaf;
    this->fullQueryStats.avgLeaf += qS.numLeaf;

    //geometric tests statistics
    if(this->fullQueryStats.minGeometricTest > qS.numGeometricTest)
        this->fullQueryStats.minGeometricTest = qS.numGeometricTest;
    if(this->fullQueryStats.maxGeometricTest < qS.numGeometricTest)
        this->fullQueryStats.maxGeometricTest = qS.numGeometricTest;
    this->fullQueryStats.avgGeometricTest += qS.numGeometricTest;

    //local analysis of access-per-tri
    int multiple = 0;
    int mult_counter = 0;
    int unique = 0;
    for(ivect_iter iter=qS.get_access_per_tri_array().begin(); iter!=qS.get_access_per_tri_array().end(); ++iter)
    {
        if(*iter==1)
            unique++;
        else
        {
            multiple += *iter;
            mult_counter++;
        }
    }

    //unique triangles access statistics
    if(this->fullQueryStats.minUniquetriAccess > unique)
        this->fullQueryStats.minUniquetriAccess = unique;
    if(this->fullQueryStats.maxUniquetriAccess < unique)
        this->fullQueryStats.maxUniquetriAccess = unique;
    this->fullQueryStats.avgUniquetriAccess += unique;

    //multiple triangles access statistics
    if(multiple != 0)
    {
        if(this->fullQueryStats.minMultipletriAccess > multiple)
            this->fullQueryStats.minMultipletriAccess = multiple;
        if(this->fullQueryStats.maxMultipletriAccess < multiple)
            this->fullQueryStats.maxMultipletriAccess = multiple;
        this->fullQueryStats.avgMultipletriAccess += multiple;
    }

    return hit_ratio;
}
