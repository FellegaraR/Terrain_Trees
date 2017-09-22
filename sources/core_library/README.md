# Terrain Trees library #

Terrain trees are a new in-core family of spatial indexes for the representation 
and analysis of Triangulated Irregular Networks (TINs).
Terrain trees combine a minimal encoding of the connectivity of the
underlying triangle mesh with a hierarchical spatial index, implicitly
representing the other topological relations among vertices, edges
and vertices. Topological relations are extracted locally within each
leaf block of the hierarchal index at runtime, based on specific application 
needs. We have developed a tool based on Terrain trees for
terrain analysis, which includes state-of-the-art estimators for slope
and curvature, and for the extraction of critical points. 
By working on TINs generated from big LiDAR (Light, Detection and Ranging) 
data sets, we demonstrate the effectiveness and scalability of the 
Terrain trees against a state-of-the-art compact data structures.

### Reference Paper ###

Riccardo Fellegara, Federico Iuricich, and Leila De Floriani. 
*Efficient representation and analysis of triangulated terrains*.
In Proceedings of SIGSPATIAL’17, Los Angeles Area, CA, USA, November 7–10, 2017, 4 pages.
[doi](https://doi.org/10.1145/3139958.3140050)

### Features ###

+ Three spatial indexes based on
    * point threshold (PR-T tree)
    * triangle threshold (PMR-T tree)
    * point and triangle thresholds (PM-T tree)
+ Two spatial decomposition
    * quadtree
    * kD-tree
+ Spatial queries
    * point location
    * box query
    * incremental nearest neighbor ([reference paper](http://link.springer.com/chapter/10.1007%2F3-540-60159-7_6))
+ Terrain Features
    * Triangle/Edges slope computation
    * Critical Points extraction
+ Curvature computation ([reference1](http://dl.acm.org/citation.cfm?id=1463498)and [reference2](http://www.umiacs.umd.edu/~deflo/papers/2010grapp/2010grapp.pdf))
    * Concentrated curvature
    * Mean and Gaussian CCurvature 
+ Soup to indexed mesh conversion
+ Points cloud indexing
    * multifield extraction

### How to compile ###

The library requires only the [boost library](http://www.boost.org/) (for dynamic_bitset class) and [cmake](https://cmake.org/) installed in your system.

Once in the root of the repository type from the command line
```
#!

cmake CMakeList.txt
```
and once configured
```
#!

make
```
This command generates a portable library file, located into *lib* folder.

The compilation process has been test on linux and mac systems.
