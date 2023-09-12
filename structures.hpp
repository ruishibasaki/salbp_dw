
#ifndef structures_hpp
#define structures_hpp

#include <vector>

struct machine{
    bool uncertain;
    bool possibly_empty;
    bool fixed_pattern;
    std::vector<int> forced;
    machine(){uncertain=false; possibly_empty=false; fixed_pattern=false;}
};


struct intpair{
    int fst, snd;
    intpair(){fst=snd=0;}
    intpair(int f,int s){fst =f; snd=s;}
    intpair& operator=(const intpair& other){fst =other.fst; snd=other.snd;}
};

struct graph_structure{
    int added_arcs;
    std::vector<int> topologic_order;
    std::vector<intpair> arcs;
    std::vector<std::vector<int> >  succ;  // Immediate successors of tasks
    std::vector<std::vector<int> >  pred;  // Immediate predecessors of tasks
    std::vector<int> pred_succ; // matrix_{nxn} (a_ij): 0 if no relation,
    //1 if j succesor of i,
    //2 if j predecessor of i.
    graph_structure(){added_arcs=0;}
    graph_structure(const graph_structure & cpy){
        topologic_order = cpy.topologic_order;
        succ = cpy.succ;
        pred = cpy.pred;
        pred_succ = cpy.pred_succ;
        arcs=cpy.arcs;
        added_arcs=0;
    }
};


#endif /* structures_hpp */
