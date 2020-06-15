#ifndef SIMPLIFICATION_AUX_STRUCTURES_H
#define SIMPLIFICATION_AUX_STRUCTURES_H

#include "basic_types/edge.h"
#include "basic_types/basic_wrappers.h"
struct Geom_Edge{
    ivect edge;
    double val;
   inline Geom_Edge(ivect &edge, double val){this->edge=edge; this->val=val;};
};

struct CompareEdge{

    bool operator()(Geom_Edge* e1,Geom_Edge* e2){

      return e1->val>e2->val;
    }
};
typedef std::priority_queue<Geom_Edge*, std::vector<Geom_Edge*>, CompareEdge> edge_queue;

class contraction_parameters
{
public: 
  contraction_parameters(){checked_edges = 0;maximum_length=0;}
    inline void increment_contracted_edges_counter() { checked_edges++; }
    inline int get_contracted_edges_num() { return checked_edges; }
    inline void set_maximum_length(double l){this->maximum_length=l;}
    inline double get_maximum_length() { return maximum_length; }
    inline void increment_counter(){simplification_counter++;}
    inline int get_counter(){return simplification_counter;}
    inline void print_simplification_counters()
    {
        if(checked_edges >0)
            cerr<<"[STATS] contracted edges "<<checked_edges<<endl;
    }
protected:
  int checked_edges;
  double maximum_length;
  int simplification_counter; //number of deleted triangles
};

#endif // SIMPLIFICATION_AUX_STRUCTURES_H