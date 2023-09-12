//
//  pattern_manager.hpp
//  
//
//  Created by Rui Shibasaki on 18/12/20.
//

#ifndef pattern_manager_hpp
#define pattern_manager_hpp

#include <stdio.h>
#include <list>
#include <vector>
#include <iostream>
#include <algorithm>
#include "data.hpp"


class Pattern{
public:
    std::vector<int> v;
    std::vector<unsigned int> image;

    bool uncertain;
    double rho;
    double cap_t;
    size_t id;
    Pattern(const std::vector<intpair>& map, const std::vector<int>& v_,
            size_t sizeOfImage, bool uncertain_, double rho_, double cap_t_, size_t id_);
    Pattern & operator = (const Pattern & pat){
        this->v=pat.v;
        this->uncertain = pat.uncertain;
        this->rho = pat.rho;
        this->cap_t = pat.cap_t;
        this->image = pat.image;
        this->id = pat.id;
        return *this;
        
    }
    
    //----------------------------------------------------------

    bool hasTask(size_t task, const std::vector<intpair>& map ) const{
        unsigned int flag = 1;
        return ( (image[map[task].fst] & (flag << (map[task].snd))) != 0 ) ;
    }
    //----------------------------------------------------------

    bool operator==(const Pattern & p1){
        for (size_t p=p1.image.size(); p--; ) {
            if (this->image[p] != p1.image[p]) {
                return false;
            }
        }
        return true;
    }

};

//===============================================================================
//===============================================================================
//===============================================================================

class PatternCollection{
public:
    std::list<Pattern> patterns;
    std::vector<int> patterns_idx; //dictionary for indices machine-pattern.
                                    //kth position indicates the index of column where its
                                    //patterns are in `patterns'.    
    
    size_t sizeOfImage;
    std::vector<intpair> map; // fst = pos, snd = bit;

    PatternCollection(size_t m, size_t n);
    size_t get_size() const{return patterns.size();}
    size_t get_qnt_k(int k){return patterns_idx[k+1]-patterns_idx[k]-1;}
};

//===============================================================================
//===============================================================================
//===============================================================================

class PatternManager{
public:

    const std::vector<int>& pred_succ; 
    const std::vector<task>& T;
    const std::vector<intpair>& qj;
    const std::vector<machine>& machines;
    
    std::vector<int> anc; //matrix for each node in line
    std::vector<int> dec; //1st col. nb of ancestors/decendents (anc/dec).
                          //2nd-(n+1)th cols place for ancestor/decendents indexes
    
    std::vector<double> certain_min_slack;
    
    const int n; const int m; const double cap;
    const int n_unc, m_unc;
    int unique_orig;
    int unique_dest;
    
    PatternCollection * collection;
    double LB;
    double UB;
    //----------------------------------------------------------
    
    PatternManager(const data &data_salbp):
        n(data_salbp.n),
        m(data_salbp.m),
        cap(data_salbp.cap),
        n_unc(data_salbp.n_unc),
        m_unc(data_salbp.m_unc),
        pred_succ(data_salbp.pred_succ),
        T(data_salbp.T),
        qj(data_salbp.qj),
        machines(data_salbp.machines)
    {
        certain_min_slack.resize(m,cap);
        UB = cap;
        LB = 0;
        unique_orig = unique_dest = -1;
    }
    
    
    //----------------------------------------------------------

    void set_lb(double lb){if(LB < lb)LB = lb;}
    void set_ub(double ub){if(UB > ub)UB = ub;}
    void set_new_ub(double ub);
    const Pattern* push_back_pattern(size_t k, bool unc, double r , const std::vector<int>& S, bool check_presence=false);
    void setCollectionPtr(PatternCollection * ptr){ collection = ptr;}
    //----------------------------------------------------------
        
    
    double localStabilityRadius(const std::vector<int>& S,
                                std::string norm,  bool machine_unc, bool & uncertain);
    
    void vector2Pattern(const std::vector<int>& solution);
    void solutions2Pattern(size_t best_id, std::vector<std::vector<int> >& solutions, std::vector<int> & best);

    bool check_valid_pattern2(const std::vector<int>& S, const std::vector<int>& forced);
    
    bool check_valid_pattern(const std::vector<int>& S,
                             const std::vector<int>& anc_dec);
 
    
    //----------------------------------------------------------
    
};

#endif /* pattern_gen_hpp */


