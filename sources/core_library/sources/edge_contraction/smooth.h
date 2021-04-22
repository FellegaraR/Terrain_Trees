#ifndef SMOOTH_H
#define SMOOTH_H

#include "terrain_trees/node_v.h"
#include "terrain_trees/node_t.h"
#include "basic_types/lru_cache.h"
#include "terrain_trees/prt_tree.h"
#include "statistics/statistics.h"
#include "utilities/container_utilities.h"
#include "terrain_trees/mesh_updater.h"
#include "utilities/usage.h"
#include "utilities/cli_parameters.h"
#include "queries/topological_queries.h"
#include "simplification_aux_structures.h"

class Smooth
{
public:
    Smooth(){}
    void smooth(PRT_Tree &tree, Mesh &mesh, double lambda, int iter_num,cli_parameters &cli);
    void smooth_parallel(PRT_Tree &tree, Mesh &mesh, double lambda, int iter_num,cli_parameters &cli);

protected:
    void smooth_compute(Node_V &n, Mesh &mesh, Spatial_Subdivision &division);

    void smooth_leaf(Node_V &n, Mesh &mesh);

    inline void set_zvalues(Mesh &mesh){
        #pragma omp parallel for
        for(int i = 1; i<=new_elevations.size();i++){
            mesh.get_vertex(i).set_c(2,new_elevations[i-1]);

        }
    }

    int iter_num;
    double lambda;
    dvect new_elevations;
    bool keep_orig=false;
};

#endif  