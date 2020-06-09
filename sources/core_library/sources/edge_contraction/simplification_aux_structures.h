#include "basic_types/edge.h"
#include "basic_types/basic_wrappers.h"
struct Geom_Edge{
    Edge* edge;
    double val;
   inline Geom_Edge(Edge* edge, double val){this->edge=edge; this->val=val;};
};

struct CompareEdge{

    bool operator()(Geom_Edge e1,Geom_Edge e2){

      return e1.val<e2.val;
    }
};
typedef std::priority_queue<Geom_Edge*, std::vector<Geom_Edge*>, CompareEdge> edge_queue;