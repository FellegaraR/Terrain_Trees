# Terrain Trees library #

We propose a new in-core family of spatial indexes for the analysis of 
Triangulated Irregular Networks (TINs). We call such indexes Terrain trees. 
A Terrain tree combines a minimal encoding of the connectivity of the 
underlying triangle mesh with a hierarchical spatial index implicitly 
encoding the other topological relations among the mesh elements. 
Topological relations are extracted locally within each leaf block of 
the hierarchy at runtime, based on specific application needs. We 
introduce a new tool for the multivariate analysis of the surface by 
combining different scalar fields (i.e., the elevation and the curvature 
values). By computing a combinatorial discrete vector field we are able 
to study the mutual relationships between the different fields providing 
new insight on the terrain morphology. Moreover, we have developed other 
state-of-the-art estimators, such as slope estimation, curvature 
computation, and the extraction of the terrain critical points on the 
Terrain trees.

### Reference Paper ###

add a reference paper here

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
This command generates a portable library file, located into *lib* folder, as well as an executable into the `bin` folder.

The compilation process has been test on linux and mac systems.

### Execute unit tests ###

To compile the unit-tests run from command line the following command
```
#!

make tests
```
Once compiled, it is possible to test the main functionalities running a script located into the `bin` folder.

From command line executing 
```
#!

sh run_tests.sh
```
checks the functionalities of the main implemented features.
The output files are saved into the `data` folder (where the input datasets are located).

It is possible to clean the output files running, from command line (from `data` folder), the following command
```
#!

sh clean_up.sh
```

### Use the main library ###

In the `bin` folder there is the main executable file named `terrain_trees` that contains the whole library. For a complete list of the command line options refer the [wiki](https://bitbucket.org/riccardo_fellegara/terrain-trees/wiki/edit/Command%20line%20parameters) page.
