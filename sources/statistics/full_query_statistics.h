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

#ifndef _FULLQUERYSTATISTICS_H
#define	_FULLQUERYSTATISTICS_H

#include <cstddef>
using namespace std;

///A class storing the statistics obtained from a series of query
class FullQueryStatistics{
public:
    ///A constructor method
    FullQueryStatistics()
    {
        mintri = 10000, maxtri=0;
        minNode = 10000, maxNode=0, minLeaf = 10000, maxLeaf=0;
        minGeometricTest = 10000, maxGeometricTest=0;
        minUniquetriAccess = 10000; maxUniquetriAccess = 0;
        minMultipletriAccess = 10000; maxMultipletriAccess = 0;
        avgtri = avgNode = avgLeaf = avgGeometricTest = avgMultipletriAccess = avgUniquetriAccess = 0;
    }

    ///A public variable representing the minimum number of triangles found during query
    utype mintri;
    ///A public variable representing the maximum number of triangles found during query
    utype maxtri;
    ///A public variable representing the minimum number of node visited during query
    int minNode;
    ///A public variable representing the maximum number of node visited during query
    int maxNode;
    ///A public variable representing the minimum number of leaf node visited during query
    int minLeaf;
    ///A public variable representing the minimum number of leaf node visited during query
    int maxLeaf;
    ///A public variable representing the minimum number of geometric test executed during query
    int minGeometricTest;
    ///A public variable representing the maximum number of geometric test executed during query
    int maxGeometricTest;
    ///A public variable representing the average number of triangles found during query
    coord_type avgtri;
    ///A public variable representing the average number of node visited during query
    coord_type avgNode;
    ///A public variable representing the average number of leaf node visited during query
    coord_type avgLeaf;
    ///A public variable representing the average number of geometric test executed during query
    coord_type avgGeometricTest;
    ///A public variable representing the minimum number of triangles that have been uniquely accessed during a query
    int minUniquetriAccess;
    ///A public variable representing the maximum number of triangles that have been uniquely accessed during a query
    int maxUniquetriAccess;
    ///A public variable representing the average number of triangles that have been uniquely accessed during a query
    coord_type avgUniquetriAccess;
    ///A public variable representing the minimum number of triangles that have been accessed several times during a query
    int minMultipletriAccess;
    ///A public variable representing the maximum number of triangles that have been accessed several times during a query
    int maxMultipletriAccess;
    ///A public variable representing the average number of triangles that have been accessed several times during a query
    coord_type avgMultipletriAccess;

};

#endif	/* _FULLQUERYSTATISTICS_H */

