//
//  pattern_model.hpp
//  
//
//  Created by Rui Shibasaki on 23/12/20.
//

#ifndef solver_hpp
#define solver_hpp

#include <stdio.h>
#include "data.hpp"
#include "pattern_manager.hpp"

#include <ilcplex/ilocplex.h>

#define ADD_VALID_INEQ 1 
#define ADD_PRECED_CONSTRAINTS 0 // 0 - do not add them; 1 - add them all; 2 - add them on the fly   
#define REINFORCE_PRICE_WITH_PRECED 1 // 0 - do not reinforce; 1 - reinforce with local precedence constraints    


struct pricing_problem{
    IloNumVarArray xis;
    IloNumVarArray x;

    #if ADD_VALID_INEQ ==1 
    IloNumVar rho;
    #endif

    IloObjective obj;
    IloModel model2;
};

class Solver{
public:
    int nfix_var;
    const int n, m, n_unc, m_unc;
    const double cap;
    double UB;
    double LB;
    const data * data_salbp;
    PatternCollection * collection;
    PatternManager * pattern_manager;
    clock_t t_u;
    double time_to_feas;

    IloEnv env;
    IloCplex cplex;
    IloCplex cplex_pricing;

    IloModel model;
    IloNumVarArray lam;
    IloNumVarArray lam0;
    
    IloNumVarArray slacks_j;
    IloNumVarArray slacks_k;
    
    IloRangeArray c1;
    IloRangeArray c2;
    IloRangeArray c3;
    IloRangeArray c4;

    #if ADD_VALID_INEQ ==1 
    IloRangeArray cuts;
    #endif

    IloNumVar rho;

    //dual solution
    IloNumArray mu;
    IloNumArray w;
    IloNumArray kappa;
    IloNumArray pi;
    IloNumArray zeta;

    std::vector<int> mapi;
    std::vector<int> fixed_vars;
   
    //subproblem
    std::vector<pricing_problem> pricing_problems;
    
    
    Solver(double lb, double ub, const data& data_salbp_, PatternManager * pattern_manager_);
   
    //================================================================

    void set_parameters();
    void create_model(const std::vector<intpair>& arcs,
                      const std::vector<intpair>&  qj,
                      const std::vector<machine>&  machines,
                      bool add_slacks);
    
    
    void solve(const std::string& instance, bool add_slacks, std::vector<int> best_id);
    
    bool find_feas_basis();
    bool solve_colgen();
    //================================================================
    void create_subproblem(const std::vector<task>& T,
                           const std::vector<machine> & machines,
                           const std::vector<int>&  tasks_in_machine);
    
    void modify_subproblems(const std::vector<intpair>& arcs,
                            const std::vector<intpair>& qj,
                            const std::vector<machine>& machines,
                            const std::vector<task>& T,
                            const std::vector<int>&  tasks_in_machine);
    
    void translate_column(size_t k, size_t position,
                          const std::vector<intpair>&qj,
                          const std::vector<intpair>& arcs,
                          const std::vector<int>&  tasks_in_machine );
    void add_col(size_t k,  IloNum rho_rk,
                 const IloNumArray & xis_,
                 const std::vector<int>& in_k,
                 const std::vector<intpair>&qj,
                 const std::vector<intpair>& arcs,
                 const Pattern * pattern );
    
    #if ADD_PRECED_CONSTRAINTS==2
    bool add_row(const std::vector<intpair>&qj,
                const std::vector<intpair>& arcs);
    #endif
    
    bool solve_pricing(double & current_ub);
    
    
    //================================================================
    void variable_fixing(double current_ub,
                         const IloNumArray &lam_,
                         const IloNumArray& rc);
    void spread_new_ub(const std::vector<machine> & machines,
                       const std::vector<int>&  tasks_in_machine);
    double solve_mip(double & time,const std::vector<intpair>&qj,  std::vector<int>& best_id);


};


#endif /* pattern_model_hpp */
