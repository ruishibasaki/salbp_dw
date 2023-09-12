//
//  heuristic.hpp
//  
//
//  Created by Rui Shibasaki on 30/04/21.
//

#ifndef heuristic_hpp
#define heuristic_hpp

#include <stdio.h>
#include "data.hpp"
#include <vector>
#include <utility>



class heuristicSolution{
public:
    
    const int n, m, n_unc, m_unc;
    const double cap;
    const double UB;
    double LB;
    data * data_salbp;
    
    size_t best_sol_id;

    std::vector<int> candidates;
    std::vector<std::pair<bool,size_t> > assigned;
    std::vector<std::pair<double, unsigned int> > mach_state;
    
    //=====================================================================

    heuristicSolution(double lb, double ub, data & data_salbp_):
        n(data_salbp_.n),
        m(data_salbp_.m),
        cap(data_salbp_.cap),
        n_unc(data_salbp_.n_unc),
        m_unc(data_salbp_.m_unc), LB(lb), UB(ub){
            data_salbp = &data_salbp_;
            candidates.reserve(n);
            assigned.resize(n,std::pair<bool,size_t>(false,0));
            mach_state.resize(m,std::pair<double, unsigned int>(cap,0));
            srand(14);
            LB = LB ==0 ? -1.0 : LB;
            best_sol_id =0;
        }
    
    //=====================================================================
    double solve(std::vector<std::vector<int> >& solutions);
    double evaluate();
    void translate_solution(std::vector<int>& solution);
    void create_candidate_list(size_t k);
    void assign_task(size_t k);

    //=====================================================================

    bool check_candidate(size_t k, size_t task);
    bool check_stability_level(size_t k, size_t task);
    bool check_assign_intvl(size_t k, size_t task);
    bool check_predecessors(size_t k, size_t task);
};



#endif /* heuristic_hpp */
