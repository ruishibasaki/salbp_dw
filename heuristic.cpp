//
//  heuristic.cpp
//  
//
//  Created by Rui Shibasaki on 30/04/21.
//

#include "heuristic.hpp"
#include <algorithm>
//=====================================================================
//=====================================================================
//=====================================================================

bool
heuristicSolution::check_stability_level(size_t k, size_t tsk){
    const std::vector<task>& T = data_salbp->T;
    const std::vector<machine>& machines = data_salbp->machines;
    double cap_t = mach_state[k].first;
    double modv = mach_state[k].second;
    cap_t -= T[tsk].time;
    if(T[tsk].is_unc() || machines[k].uncertain)
        modv+= 1.0;
    
    if (modv>0) {
        cap_t /= modv;
    }
    if(cap_t>=LB && cap_t>=0 ){
        return true;
    }else return false;
}

//=====================================================================

bool
heuristicSolution::check_assign_intvl(size_t k, size_t task){
    const std::vector<intpair>& qj = data_salbp->qj;
    ++k;
    if(qj[task].fst<=k && k<=qj[task].snd){
        return true;
    }else return false;
}

//=====================================================================


bool
heuristicSolution::check_predecessors(size_t k, size_t task){
    for (size_t j=n; j--; ) {
        if ((data_salbp->pred_succ[task*n+j]==2) && !assigned[j].first) {
            return false;
        }
    }
    return true;
}

//=====================================================================

bool
heuristicSolution::check_candidate(size_t k, size_t j){
    if(! check_assign_intvl(k,j)){
        return false;
    }
    if (!check_stability_level(k,j)) {
        return  false;
    }
    if (!check_predecessors(k,j)) {
        return false;
    }
    return true;
}

//=====================================================================
//=====================================================================
//=====================================================================

void
heuristicSolution::create_candidate_list(size_t k){
    candidates.clear();
    for (size_t j=n; j--; ) {
        if (!assigned[j].first && check_candidate(k,j)) {
            candidates.push_back(j);
        }
    }
}

//=====================================================================

void
heuristicSolution::assign_task(size_t k){
    const std::vector<task>& T = data_salbp->T;
    const std::vector<machine>& machines = data_salbp->machines;
    
    std::random_shuffle(candidates.begin(), candidates.end());
    size_t tsk = rand()%candidates.size();
    tsk = candidates[tsk];
    assigned[tsk].first = true;
    assigned[tsk].second = k;
    mach_state[k].first -= T[tsk].time;
    if(T[tsk].is_unc() || machines[k].uncertain)
        ++mach_state[k].second;
    
    //std::cout<<"task "<<tsk+1<<" machine "<<k+1<<" rho: "<<mach_state[k].first<<" / "<<mach_state[k].second<<std::endl;
}

//=====================================================================

double
heuristicSolution::solve(std::vector<std::vector<int> >& solutions){
    const std::vector<task>& T = data_salbp->T;
    const std::vector<machine>& machines = data_salbp->machines;
    const std::vector<intpair>& qj = data_salbp->qj;

    size_t max_iters = 5;
    size_t it =0;
    size_t sol_n =0;
    double rho_it ;
    double rho=-1;
    size_t maxsol= solutions.size();
    while (it<5000) {
        for (size_t j=n; j--; ) {
            if(qj[j].fst==qj[j].snd){
                size_t k = qj[j].fst-1;
                assigned[j].first = true;
                assigned[j].second = k;
                mach_state[k].first -= T[j].time;
                if(T[j].is_unc() || machines[k].uncertain)
                    ++mach_state[k].second;
            }
        }
        
        for (size_t k =0; k<m; ++k) {
            while(1){
                create_candidate_list(k);
                if(candidates.empty()){
                    break;
                }else assign_task(k);
            }
        }
        rho_it = evaluate();
        if(rho_it>=0 && rho_it>=rho){
            //std::cout<<it<<" rho_it: "<<rho_it<<std::endl;
            if(rho_it>rho){
                it=0;
                rho = rho_it;
                LB = rho_it;
                best_sol_id = sol_n;
                data_salbp->reprocess(rho);
            }
            translate_solution(solutions[sol_n++]);
            if (sol_n>=maxsol){sol_n=0;}
        }
        assigned.assign(n,std::pair<bool,size_t>(false,0));
        mach_state.assign(m,std::pair<double, unsigned int>(cap,0));
        ++it;
    }
    
    return rho;
}

//=====================================================================


double
heuristicSolution::evaluate(){
    size_t k;
    for (size_t j=n; j--; ) {
        if (!assigned[j].first) {
            //std::cout<<j+1<<" not assigned "<<std::endl;
            return -1;
        }
    }
    double rho = cap;
    for (k =0; k<m; ++k) {
        if(mach_state[k].second==0) continue;
        double cap_t = mach_state[k].first;
        double modv =  mach_state[k].second;
        cap_t /=  modv;
        //std::cout<<"machine "<<k+1<<" rho: "<<cap_t<<std::endl;
        rho  = rho > cap_t ? cap_t : rho;
    }
    return rho;
}

//=====================================================================

void
heuristicSolution::translate_solution(std::vector<int>& solution){
    solution.assign(n*m,0);
    size_t k;

    for (size_t j=n; j--; ) {
        if (assigned[j].first) {
            k = assigned[j].second;
            solution[k*n+j]=1;
            //std::cout<<"("<<j+1<<","<<k+1<<") ";
        }else{
            //std::cout<<j+1<<" not assigned "<<std::endl;
            return;
        }
    }
    //std::cout<<std::endl;
}

//=====================================================================
//=====================================================================
//=====================================================================
