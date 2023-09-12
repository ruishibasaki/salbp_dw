
#ifndef data_hpp
#define data_hpp

#include "task.hpp"
#include "structures.hpp"
#include <vector>
#include <list>
#include <math.h>
#include <stdlib.h>
#include <algorithm>    // std::sort
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

class data{
public:
    int   n;                     // The number of tasks
    int   m;                     // The number of bins
    double cap;                   // The capacity of bin
    
    int   n_unc;                 // The number of uncertain tasks
    int   m_unc;                 // The number of uncertain workstations
    
    int unique_orig;
    int unique_dest;
    
    std::vector<task> T;              // The tasks
    std::vector<machine> machines;       // The machines - uncertain or not
    

    graph_structure graph_struct;
    //----------------alias--------------
    std::vector<int>& topologic_order = graph_struct.topologic_order;
    std::vector<std::vector<int> >&  succ = graph_struct.succ;  // Immediate successors of tasks
    std::vector<std::vector<int> >&  pred = graph_struct.pred;  // Immediate predecessors of tasks
    std::vector<int>& pred_succ = graph_struct.pred_succ; // matrix_{nxn} (a_ij): 0 if no relation,
    //1 if j succesor of i,
    //2 if j predecessor of i.
    std::vector<intpair>& arcs = graph_struct.arcs;  // Arcs
    //------------------------------

    std::vector<intpair>qj;  // Intervals Q(j) for each task j.
    std::vector<int>tasks_in_machine; //n_k: the number of tasks in machine k
                                     // position tasks_in_machine[n_k] = n_{k+1}
                                     // tasks_in_machine[0] = n_0
    
    std::vector<int> perm_task;       // Permutation of tasks
    std::vector<int> perm_bin;        // Permutation of bins
    
    data(){
        n=m=n_unc=m_unc=0;
        cap=0;
        unique_orig =unique_dest=-1;
    }
   
    //==========================================================
    
    void read_data(const char* const argv[]);
    bool test_sols_validity(size_t& best_id, std::vector<std::vector<int> >& solutions);

    void compute_unique_orig_dest();
    void make_sets_taks_machine();
    void compute_topologic_order(graph_structure& gs);
    
    bool compute_intervals(double LB, std::vector<intpair>& qj_temp, graph_structure& gs);
    bool compute_intervals_PatersonAlbracht(double LB, std::vector<intpair>& qj_temp) const;
    bool compute_intervals_Pirogov(double LB, std::vector<intpair>& qj_temp) const;
    bool compute_intervals_1rCmax(double LB, std::vector<intpair>& qj_temp, graph_structure& gs);
    double find_makespan(double LB, const std::vector<double> & est, 
                    std::vector<std::pair<double, int> > & tasks) const;
    double improve_est(double t, double est, double LB, const std::vector<int>& pred) const;
    void add_new_arc(int from, int to, graph_structure& gs);


    void preprocess(double LB);
    void reprocess(double LB);
    
    double compute_ub( double UB);
    double improve_ub( double UB, double LB);
};




#endif /* structures_hpp */
