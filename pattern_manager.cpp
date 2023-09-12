//
//  pattern_manager.cpp
//  
//
//  Created by Rui Shibasaki on 18/12/20.
//

#include "pattern_manager.hpp"
#include <list>
#include <cmath>



//============================================================
//============================================================
//============================================================

Pattern::Pattern(const std::vector<intpair>& map, const std::vector<int>& v_,
                 size_t sizeOfImage, bool uncertain_, double rho_, double cap_t_, size_t id_){
    v = v_;
    uncertain = uncertain_;
    rho = rho_;
    cap_t = cap_t_;
    id = id_;
    image.resize(sizeOfImage); 
    unsigned int flag = 1;
    size_t task;
    size_t sz = v.size();
    for(;sz--;){
        task = v[sz]-1;
        // Set the bit at the k-th position in image[map[arc].fst]
        image[map[task].fst] |= (flag << map[task].snd);
    }
}

//============================================================
//============================================================
//============================================================

PatternCollection::PatternCollection(size_t m, size_t n){
    patterns_idx.resize(m+1,0);
    size_t sizeOfInt=8*sizeof(unsigned int);
    sizeOfImage = (n/sizeOfInt)+1;

    map.resize(n);
    for(size_t i=0;i<n;++i){
        map[i].fst = i/sizeOfInt;
        map[i].snd = i%sizeOfInt;
        //std::cout<<i<<" "<<map[i].fst<<" : "<<map[i].snd<<std::endl;
    }
}

//================================================================
//================================================================
//================================================================

const Pattern*
PatternManager::push_back_pattern(size_t k, bool unc, double r , const std::vector<int>& S, bool check_presence ){
    size_t begin = collection->patterns_idx[k];
    std::list<Pattern>::iterator it = collection->patterns.begin();
    std::advance(it,begin);
    Pattern * pattern = new Pattern(collection->map, S, collection->sizeOfImage, unc,  (unc ? std::min(r,UB): UB), r, collection->get_size() );
    //std::cout<<"rho :"<<pattern->rho<<std::endl;
    if (check_presence) {
        std::list<Pattern>::const_iterator it2 = collection->patterns.begin();
        std::advance(it2,begin);
        for (size_t l = collection->patterns_idx[k]; l<collection->patterns_idx[k+1] ; ++l){
            if(*pattern == (*it2)){
                delete pattern;
                return &(*it2);
            }
            ++it2;
        }
    }
    
    collection->patterns.insert(it, *pattern );
    
    for (size_t kk=k; kk<m; ++kk) {
        ++collection->patterns_idx[kk+1];
    }
    return pattern;
}

//================================================================

void
PatternManager::set_new_ub(double ub){
    set_ub(ub);
    size_t szcollect = collection->get_size();
    std::list<Pattern>::iterator it = collection->patterns.begin();
    for (size_t l = 0; l<szcollect ; ++l){
        if((*it).rho > UB){(*it).rho = UB;}
        ++it;
    }
}

//================================================================
//================================================================
//================================================================


//returns the local stability radius of a pattern S, in norm norm.
double
PatternManager::localStabilityRadius(const std::vector<int>& S,
                                 std::string norm, bool machine_unc,  bool & uncertain){
    
    int ts = 0;
    int n = 0;
    int i;
    for(int p =S.size(); p--;){
        i = S[p]-1;
        ts += T[i].time;
        if(T[i].is_unc()) ++n;
    }
    if(machine_unc) n = S.size();
    else if(n==0){n=1;}
    else uncertain=true;
    
    if(norm == "Linf")
        return (cap - ts) /double(n);
    else
        return cap - ts;
}


//================================================================
//Expect a solution ordered by machines {0..m}, 
//-1 after the last task assigned to the machine
void 
PatternManager::vector2Pattern(const std::vector<int>& solution){
    bool machine_unc, uncertain;
    double rho_local;
    int id=0;
    int task;
    std::vector<int> S;
    S.reserve(n);
    for (size_t k = 0; k < m; ++k){
        machine_unc = machines[k].uncertain;
        uncertain = machine_unc;
        task = solution[id++];
        while(task>0){
            S.push_back(task);
            task = solution[id++];
        }
        if(S.size()==0) continue;

        rho_local = localStabilityRadius( S, "Linf", machine_unc,  uncertain);
        push_back_pattern(k, uncertain, rho_local, S );
        S.clear();
    }
}

//================================================================
//Expect a solution {m*n} solution[k*n+i]=1 if taks i is assigned to machine k,
//0 otherwise
void
PatternManager::solutions2Pattern(size_t best_id, std::vector<std::vector<int> >& solutions, std::vector<int> & best){
    bool machine_unc, uncertain, equal;
    double rho_local;
    int task;
    for (size_t sz=solutions.size(); sz--; ) {
        if (solutions[sz].empty()){ continue;}

        std::vector<int> S;
        S.reserve(n);
        //std::cout<<"solution: "<<std::endl;
        for (size_t k = 0; k < m; ++k){
            //std::cout<<"station "<<k+1<<std::endl;
            machine_unc = machines[k].uncertain;
            uncertain = machine_unc;
            for (size_t j = 0; j < n; ++j){
                if (solutions[sz][k*n+j]==1) {
                    //std::cout<<j+1<<" ";
                    S.push_back(j+1);
                }
            }
            if(S.size()==0) continue;
            
            rho_local = localStabilityRadius( S, "Linf", machine_unc,  uncertain);
            const Pattern* pattern = push_back_pattern(k, uncertain, rho_local, S, true );
            if(pattern!=0 && best_id==sz)
                best.push_back(pattern->id);
            S.clear();
            //std::cout<<" rho: "<<rho_local<<std::endl;
        }
    }
}

//================================================================
//Check validity for patterns in machines 1 and m. Check if all precedents(descendants) of jobs
//in those machine patterns are present in the pattern.
bool
PatternManager::check_valid_pattern(const std::vector<int>& S,
                                  const std::vector<int>& anc_dec){
    
    int szanc_dec;
    int sz_s = S.size();
    int i, idx;
    bool ensured=true;
    for(int p=0;p<sz_s;++p){
        i = S[p];
        idx = (i-1)*(n+1);
        szanc_dec = anc_dec[idx];
        ++idx;
        if(szanc_dec==0) continue;
        for(int l = 0;l<szanc_dec;++l){
            int job = anc_dec[idx++];
            ensured=false;
            for(int j=0;j<sz_s;++j){
                if(job == S[j]){ensured=true;  break;}
            }
            if(!ensured) return false;
        }
    }
    return true;
    
}

//================================================================
//Check validity for patterns in machines 2 to m-1. Check if machine-forced jobs are present
//in all patterns for that machine.

bool
PatternManager::check_valid_pattern2(const std::vector<int>& S, const std::vector<int>& forced){
    if(forced.empty()) return true;
    int sz_s = S.size();
    int sz_f = forced.size();
    int i;
    bool ensured;
    for (int f=0; f<sz_f; ++f) {
        i = forced[f];
        ensured=false;

        for(int p=0;p<sz_s;++p){
            if(i == S[p]) ensured = true;
        }
        if(!ensured) return false;
    }
    
    return true;
}
