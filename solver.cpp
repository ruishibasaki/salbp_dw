//
//  pattern_model.cpp
//  
//
//  Created by Rui Shibasaki on 23/12/20.
//


#include "solver.hpp"
#include <time.h>

#define TIME_LIMIT 10.0
#define EPSILON 5e-4
#define PRES 1e-10
#define M 1e10


//================================================================
//================================================================
//================================================================

ILOMIPINFOCALLBACK0(InfoMIPCallBack){
    
    if(getCplexTime() - getStartTime() > 10.0){
        if (getBestObjValue() < 0) {
            abort();
            return;
        }else if(hasIncumbent() && (getIncumbentObjValue()) > EPSILON){
            abort();
            return;
        }else if(getCplexTime() - getStartTime() > TIME_LIMIT){
            abort();
            return;
        }
    }
    
}

//================================================================

bool verify_timeout(bool & time_out, clock_t t_u ){
    double t = double( clock() - t_u ) / double( CLOCKS_PER_SEC );
    if(t > TIME_LIMIT){
        time_out = true;
        return true;
    }
    return false;
}
//================================================================
//================================================================
//================================================================

Solver::Solver(double lb, double ub, const data& data_salbp_, PatternManager * pattern_manager_):
LB(lb),
UB(ub),
n(data_salbp_.n),
m(data_salbp_.m),
cap(data_salbp_.cap),
n_unc(data_salbp_.n_unc),
m_unc(data_salbp_.m_unc),
cplex(env), model(env), cplex_pricing(env),
c1(env),c2(env), c3(env), c4(env), 
mu(env), w(env), kappa(env), pi(env)
#if ADD_VALID_INEQ ==1 
,cuts(env), zeta(env)
#endif
{
    data_salbp =&data_salbp_;
    pattern_manager = pattern_manager_;
    collection = pattern_manager_->collection;
    pricing_problems.resize(m);
    mapi.resize(data_salbp_.arcs.size()*m, -1);
    time_to_feas=0.0;
    srand(111);
}

//================================================================
//================================================================
//================================================================
void
Solver::set_parameters(){
    cplex.setParam(IloCplex::Param::Threads,1);
    cplex.setParam(IloCplex::Param::RootAlgorithm, 1);
    //cplex.setParam(IloCplex::PreInd, 0);
    cplex.setParam(IloCplex::Param::ClockType, 1);
    //cplex.setParam(IloCplex::Param::MIP::Limits::EachCutLimit, 0);
    //cplex.setParam(IloCplex::Param::MIP::Cuts::Gomory, -1);
    //cplex.setParam(IloCplex::Param::MIP::Cuts::LiftProj, -1);
    //cplex.setParam(IloCplex::MIPDisplay, 4);
    //cplex.setParam(IloCplex::NodeLim, 10);
    cplex.setOut(env.getNullStream());
    //cplex.setParam(IloCplex::Param::TimeLimit, 3600.0); // Time limit in seconds
    cplex.setParam(IloCplex::Param::Advance, 2);
    
    //cplex_pricing.setParam(IloCplex::Param::Advance, 0);
    //cplex_pricing.setParam(IloCplex::Param::TimeLimit, 10.0); // Time limit in seconds
    cplex_pricing.setParam(IloCplex::Param::Threads,1);
    cplex_pricing.setOut(env.getNullStream());
    cplex_pricing.setParam(IloCplex::Param::ClockType, 1);
    //cplex_pricing.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, PRES);
    cplex_pricing.use(InfoMIPCallBack(env));
}

//================================================================

void
Solver::create_model(const std::vector<intpair>& arcs,
                     const std::vector<intpair>&  qj,
                     const std::vector<machine>&  machines,bool add_slacks){
    
    int szarcs = arcs.size();
    size_t szcollect = collection->get_size();
    lam = IloNumVarArray(env, szcollect, 0.0, IloInfinity, ILOFLOAT); //
    lam0 = IloNumVarArray(env, m, 0.0, IloInfinity, ILOFLOAT); //
    
    rho = IloNumVar(env);
    
    if(add_slacks){
        slacks_j = IloNumVarArray(env, n, 0.0, IloInfinity, ILOFLOAT);
        slacks_k = IloNumVarArray(env, m, 0.0, IloInfinity, ILOFLOAT);
    }
    
    IloExpr obj(env);
    obj += rho;
    if(add_slacks){
        for(size_t k =0; k<m; ++k) {
            obj += -1e4*slacks_k[k];
        }
        for(size_t j =0; j<n; ++j) {
            obj += -1e4*slacks_j[j];
        }
    }
    
    int cont_lam0=0;
    int cont=0;
    

    std::list<Pattern>::const_iterator it;
    for(size_t k =0; k<m; ++k) {
        IloExpr constraint(env);
        it = collection->patterns.begin();
        std::advance(it,collection->patterns_idx[k]);
        for (size_t l = collection->patterns_idx[k]; l<collection->patterns_idx[k+1] ; ++l){
            constraint += lam[it->id];
            ++it;
        }
        if(add_slacks) constraint += slacks_k[k];
        constraint -= 1;
        
        if(machines[k].uncertain || 
            (collection->patterns_idx[k+1]==collection->patterns_idx[k])){
            //std::cout<<"put lam0 for machine "<<k+1<<std::endl;
            constraint += lam0[cont_lam0];
            ++cont_lam0;
        }
        c1.add((constraint == 0));
        model.add(c1[cont]);
        constraint.end();

        #if ADD_VALID_INEQ ==1 
        IloExpr constraint2(env);
        constraint2+= rho - UB;
        it = collection->patterns.begin();
        std::advance(it,collection->patterns_idx[k]);
        for (size_t l = collection->patterns_idx[k]; l<collection->patterns_idx[k+1] ; ++l){
            if(it->uncertain) constraint2 -= (it->rho - UB)*lam[it->id];
            ++it;
        }
        cuts.add((constraint2 <= 0));
        model.add(cuts[cont]);
        constraint2.end();
        #endif 

        ++cont;
   
       
    }

    cont=0;
    for(size_t j =0; j<n; ++j) {
        IloExpr constraint(env);
        it = collection->patterns.begin();
        for (size_t l = 0; l<szcollect ; ++l){
            if(it->hasTask(j, collection->map)){constraint += lam[it->id];}
            ++it;
        }
        if(add_slacks) constraint += slacks_j[j];
        constraint -= 1;
        c2.add((constraint == 0));
        model.add(c2[cont]); //model.add(constraint == 0) for sum a*lambda = 1;
        constraint.end();
        
        IloExpr constraint2(env);
        constraint2 += rho;
        for(size_t k = qj[j].fst; k<=qj[j].snd; ++k){
            it = collection->patterns.begin();
            std::advance(it,collection->patterns_idx[k-1]);
            for (size_t l = collection->patterns_idx[k-1]; l<collection->patterns_idx[k] ; ++l){
                if(it->hasTask(j, collection->map)){
                    constraint2 -= (it->rho)*lam[it->id];
                }
                ++it;
            }
        }
        c3.add((constraint2 <= 0));
        model.add(c3[cont]);
        constraint2.end();
        ++cont;

    }
    //Precedence constraints.
    #if ADD_PRECED_CONSTRAINTS ==1
    int i,j;
    cont=0;
    for(size_t a =0; a<szarcs;++a){
        i = arcs[a].fst-1;
        j = arcs[a].snd-1;
        //std::cout<<"arc "<<a+1<<" ("<<i+1<<","<<j+1<<")"<<" intv ["<<qj[j].fst<<","<<qj[i].snd<<"] nbcons "<<qj[i].snd-qj[j].fst+1<<std::endl;
        for(size_t k = qj[j].fst; k<=qj[i].snd; ++k){
            IloExpr constraint(env);
            //std::cout<<"station "<<k<<std::endl;
            for(size_t l = k-1; l<qj[i].snd; ++l){
                it = collection->patterns.begin();
                std::advance(it,collection->patterns_idx[l]);
                for (int p = collection->patterns_idx[l]; p<collection->patterns_idx[l+1] ; ++p){
                    if(it->hasTask(i, collection->map)){constraint += lam[it->id]; }
                    ++it;
                }
            }
            for(size_t l = k-1; l<qj[j].snd; ++l){
                it = collection->patterns.begin();
                std::advance(it,collection->patterns_idx[l]);
                for (int p = collection->patterns_idx[l]; p<collection->patterns_idx[l+1] ; ++p){
                    if(it->hasTask(j, collection->map)){constraint -= lam[it->id];  }
                    ++it;
                }
            }
            
           
            c4.add((constraint <= 0));
            model.add(c4[cont]);
            mapi[a*m+k-1] = cont;
            constraint.end();
            
            ++cont;
        }
    }
    #endif
    model.add(IloMaximize(env, obj));
    obj.end();
}

//================================================================
//================================================================
//================================================================

void
Solver::create_subproblem(const std::vector<task>& T,
                          const std::vector<machine> & machines,
                       const std::vector<int>&  tasks_in_machine){
    
    size_t position=0;
    size_t task;
    int number_taks_machine;
    for (size_t k=0; k<m; ++k) {
        pricing_problem& p =pricing_problems[k];
        number_taks_machine = tasks_in_machine[position];
        
        p.xis = IloNumVarArray(env, number_taks_machine, 0.0, IloInfinity, ILOFLOAT );
        p.x = IloNumVarArray(env, number_taks_machine, 0.0, 1.0, ILOINT);
        p.model2 = IloModel(env);
        
        IloExpr objexpr(env);
        for (size_t j=0; j<number_taks_machine; ++j) {        
            objexpr += p.x[j] + p.xis[j];
        }

        #if ADD_VALID_INEQ ==1 
        p.rho = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
        objexpr += p.rho;
        #endif
        
        p.obj = IloMaximize(env, objexpr);
        p.model2.add(p.obj);
        objexpr.end();
        
        IloExpr constraint(env);
        for (size_t j=0; j<number_taks_machine; ++j) {
            IloExpr constraint2(env);
            constraint2 += p.xis[j] - UB*p.x[j];
            p.model2.add((constraint2 <= 0));
            constraint2.end();
        
            #if ADD_VALID_INEQ ==1 
            IloExpr constraint3(env);
            constraint3 += p.xis[j] - p.rho;
            p.model2.add((constraint3 <= 0));
            constraint3.end();
            #endif

            task = tasks_in_machine[position+j+1]-1;
            constraint += T[task].time * p.x[j];
            if(data_salbp->qj[task].fst == data_salbp->qj[task].snd){
                p.x[j].setLB(1.0);
            }

            if (machines[k].uncertain || T[task].is_unc()) {
                constraint += p.xis[j];
            
                #if ADD_VALID_INEQ ==1 
                IloExpr constraint4(env);
                constraint4 += p.xis[j] - p.rho + UB*(1 - p.x[j]);
                p.model2.add((constraint4 >= 0));
                constraint4.end();
                #endif

            }
        }
        constraint -= cap;
        p.model2.add((constraint <= 0));
        constraint.end();

        #if REINFORCE_PRICE_WITH_PRECED==1
        size_t i,j;
        size_t szarcs = data_salbp->arcs.size();
        for(size_t a =0; a<szarcs;++a){
            i = data_salbp->arcs[a].fst-1;
            j = data_salbp->arcs[a].snd-1;
            if(data_salbp->qj[i].fst-1 > k || data_salbp->qj[i].snd-1 < k ) continue;
            if(data_salbp->qj[j].fst-1 > k || data_salbp->qj[j].snd-1 < k ) continue;
            size_t posi,posj;
            for (size_t f=0; f<number_taks_machine; ++f) {
                task = tasks_in_machine[position+f+1]-1;
                if(task == i) posi = f;
                else if(task == j) posj=f;
            }
            for (size_t f=0; f<number_taks_machine; ++f) {
                task = tasks_in_machine[position+f+1]-1;
                if(data_salbp->pred_succ[j*n+task] != 1) continue;
                IloExpr constraint5(env);
                constraint5 += p.x[posi] + p.x[f] -1.0 - p.x[posj];
                p.model2.add((constraint5 <= 0));
                constraint5.end();
            }  
        } 
        #endif
        position += number_taks_machine+1;
    }
}

//================================================================

void
Solver::modify_subproblems( const std::vector<intpair>& arcs,
                           const std::vector<intpair>& qj,
                           const std::vector<machine>& machines,
                           const std::vector<task>& T,
                           const std::vector<int>&  tasks_in_machine){
    size_t position=0;
    size_t narcs = arcs.size();
    int number_taks_machine;
    int id;
    for (size_t k=0; k<m; ++k) {
        pricing_problem& p =pricing_problems[k];
        number_taks_machine = tasks_in_machine[position];
        
        IloExpr objexpr(env); 
        objexpr -= mu[k]; //std::cout<<"cnst "<<-(mu[k] + UB*zeta[k])<<std::endl;
        #if ADD_VALID_INEQ ==1 
        objexpr += zeta[k]*p.rho - UB*zeta[k];
        #endif

        std::vector<int>in_k(n, -1);
        for (size_t j=0; j<number_taks_machine; ++j) {
            id = tasks_in_machine[position+j+1]-1;
            in_k[id]=j;
            objexpr -= (kappa[id])*p.x[j];
            objexpr += w[id]*p.xis[j];
            //std::cout<<"w[]"<<w[id]<<" "<<delta3[id]<<std::endl;
            //std::cout<<"kappa[]"<<kappa[id]<<" "<<delta2[id]<<std::endl;
        }
        
        #if ADD_PRECED_CONSTRAINTS > 0
        for (size_t a=0; a<narcs; ++a) {
            size_t i = arcs[a].fst-1;
            size_t j = arcs[a].snd-1;
            size_t b = ((qj[i].snd-1) < k) ? (qj[i].snd-1) : k ;
            //std::cout<<"arc "<<i+1<<" "<<j+1<<" qi: "<<(qj[i].fst)<<","<<qj[i].snd<<" qj: "<<(qj[j].fst)<<","<<qj[j].snd<<" k: "<<k+1<<std::endl;
            if((qj[j].fst-1) <= k && k < qj[i].snd){
                for(size_t kk = qj[j].fst-1; kk<=k; ++kk){
                    //std::cout<<kk+1<<std::endl;
                    id = mapi[a*m+kk];
                    if(id>=0){
                        int id2 = in_k[i];
                        if(id2>=0){objexpr -= pi[id]*p.x[id2];
                            //std::cout<<"pi[]"<<pi[id]<<" "<<delta4[id]<<std::endl;
                        }
                    }
                }
            }
            //std::cout<<"side"<<std::endl;
            if((qj[j].fst-1) <= k && k < qj[j].snd){
                for(size_t kk = qj[j].fst-1; kk<=b; ++kk){
                    //std::cout<<kk+1<<std::endl;
                    id = mapi[a*m+kk];
                    if(id>=0){
                        int id2 = in_k[j];
                        if(id2>=0){objexpr += pi[id]*p.x[id2];
                            //std::cout<<"pi[]"<<pi[id]<<" "<<delta4[id]<<std::endl;
                        }
                    }
                }
            }
        }
        #endif
        p.obj.setExpr(objexpr);
        objexpr.end();
        position += number_taks_machine+1;
    }
}

//================================================================

void
Solver::translate_column(size_t k, size_t position,
                         const std::vector<intpair>&qj,
                         const std::vector<intpair>& arcs,
                         const std::vector<int>&  tasks_in_machine ){
    
    size_t id;
    size_t number_taks_machine = tasks_in_machine[position];
    const pricing_problem &p = pricing_problems[k];

    IloNumArray x_(env);
    IloNumArray xis_(env);

    std::vector<int> S;
    double rho_=0;;
    int n_unc_tasks = 0;
    bool unc=data_salbp->machines[k].uncertain;
    cplex_pricing.getValues(p.x, x_);
    cplex_pricing.getValues(p.xis, xis_);

    #if ADD_VALID_INEQ ==1 
    IloNum rho_rk = cplex_pricing.isExtracted(p.rho)? cplex_pricing.getValue(p.rho):UB;
    #else
    IloNum rho_rk = UB;
    #endif

    rho_rk = rho_rk <0 ? 0.0: rho_rk;
    //std::cout<<"pattern: ";
    std::vector<int>in_k(n, -1);
    for (size_t j=0; j<number_taks_machine; ++j) {
        id = tasks_in_machine[position+j+1]-1;
        in_k[id]=j;
        if(x_[j]>0.49){
            rho_ += data_salbp->T[id].time;
            if (data_salbp->T[id].is_unc()) {
                unc = true;
                ++n_unc_tasks;
            }
            S.push_back(id+1);
            //std::cout<<xis_[j]<<" ";
            //std::cout<<id+1<<" ";
        }
    }
    //std::cout<<std::endl;
    if(data_salbp->machines[k].uncertain) n_unc_tasks = S.size();
    if(n_unc_tasks==0){n_unc_tasks=1;}
    rho_ = (cap - rho_)/double(n_unc_tasks);
    
    const Pattern * pattern = pattern_manager->push_back_pattern(k, unc, rho_ , S);
    //if(pattern->uncertain)std::cout<<"new pattern machine: "<<k+1<<" rho: "<<pattern->rho<<" "<<rho_rk<<std::endl;

    add_col(k, rho_rk, xis_,in_k, qj,arcs, pattern );
    x_.end();
    xis_.end();
}

//================================================================

void
Solver::add_col(size_t k, IloNum rho_rk,
                const IloNumArray & xis_,
                const std::vector<int>& in_k,
                const std::vector<intpair>&qj,
                const std::vector<intpair>& arcs,
                const Pattern * pattern ){
    IloNumColumn col(env);
    for (size_t kk=0; kk<m; ++kk) {
    
        if (kk==k) {
            col += c1[kk](1.0);
            #if ADD_VALID_INEQ ==1 
            if(pattern->uncertain)
                col += cuts[kk](-(rho_rk - UB));
            else col += cuts[kk](0.0);
            #endif
        }else{
            col += c1[kk](0.0);
            #if ADD_VALID_INEQ ==1 
            col += cuts[kk](0.0);
            #endif 
        }
    

    }
    for(size_t j =0; j<n; ++j) {
        if(pattern->hasTask(j, collection->map)){
            col+=c2[j](1.0);
            col+=c3[j](-rho_rk);
        }
        else{
            col+=c2[j](0.0);
            col+=c3[j](0.0);
        }
    }
    
    #if ADD_PRECED_CONSTRAINTS > 0
    int id=0;
    bool i_is_in, j_is_in;
    double mult=0;
    size_t szarcs = arcs.size();
    size_t i,j;
    for(size_t a =0; a<szarcs;++a){
        i = arcs[a].fst-1;
        j = arcs[a].snd-1;
        i_is_in = pattern->hasTask(i, collection->map);
        j_is_in = pattern->hasTask(j, collection->map);
        for(size_t kk = qj[j].fst; kk<=qj[i].snd; ++kk){
            mult=0;
            if(k >= (kk-1) && k < qj[i].snd && i_is_in) {
                mult += 1.0;
            }
            if(k >= (kk-1) && k < qj[j].snd && j_is_in) {
                mult -= 1.0;
            }
            id = mapi[a*m+kk-1];
            if(id>=0)col += c4[id](mult);
        }
    }
    #endif
    lam.add(IloNumVar(col, 0.0, IloInfinity, ILOFLOAT));
    col.end();
}

//================================================================


#if ADD_PRECED_CONSTRAINTS==2
bool
Solver::add_row(const std::vector<intpair>&qj,
                const std::vector<intpair>& arcs){

    int i,j;
    bool arc_added =false;
    size_t szarcs = arcs.size();
    size_t cont=c4.getSize();
    std::list<Pattern>::const_iterator it;

    IloNumArray lam_(env);
    cplex.getValues(lam,lam_);
    for(size_t a =0; a<szarcs;++a){
        i = arcs[a].fst-1;
        j = arcs[a].snd-1;
        //std::cout<<"arc "<<a+1<<" ("<<i+1<<","<<j+1<<")"<<" intv ["<<qj[j].fst<<","<<qj[i].snd<<"] nbcons "<<qj[i].snd-qj[j].fst+1<<std::endl;
        for(size_t k = qj[j].fst; k<=qj[i].snd; ++k){
            if(mapi[a*m+k-1]>=0) continue;
            IloExpr constraint(env);
            double viol=0.0;
            //std::cout<<"station "<<k<<std::endl;
            for(size_t l = k-1; l<qj[i].snd; ++l){
                it = collection->patterns.begin();
                std::advance(it,collection->patterns_idx[l]);
                for (int p = collection->patterns_idx[l]; p<collection->patterns_idx[l+1] ; ++p){
                    if(it->hasTask(i, collection->map)){
                        constraint += lam[it->id]; 
                        viol += lam_[it->id]; 
                    }
                    ++it;
                }
            }
            for(size_t l = k-1; l<qj[j].snd; ++l){
                it = collection->patterns.begin();
                std::advance(it,collection->patterns_idx[l]);
                for (int p = collection->patterns_idx[l]; p<collection->patterns_idx[l+1] ; ++p){
                    if(it->hasTask(j, collection->map)){
                        constraint -= lam[it->id]; 
                        viol -= lam_[it->id]; 
                    }
                    ++it;
                }
            }
            
            if(viol >= PRES){
                c4.add((constraint <= 0));
                model.add(c4[cont]);
                mapi[a*m+k-1] = cont;
                arc_added=true;
                ++cont;
            }
            constraint.end(); 
        }
    }
    lam_.end();
    return arc_added;
}
#endif

//================================================================
//================================================================
//================================================================

bool
Solver::solve_pricing(double & current_ub){
    cplex.getDuals(mu,c1);
    cplex.getDuals(kappa,c2);
    cplex.getDuals(w,c3);
    #if ADD_PRECED_CONSTRAINTS > 0
    cplex.getDuals(pi,c4);
    #endif
    #if ADD_VALID_INEQ ==1 
    cplex.getDuals(zeta,cuts);
    #endif

    //std::cout<<" rho_inf: "<<current_ub;//<<std::endl;
    //std::cout<<"modify_subproblems.."<<std::endl;
    modify_subproblems(data_salbp->arcs, data_salbp->qj,data_salbp->machines, data_salbp->T, data_salbp->tasks_in_machine);
    
    size_t position=0;
    size_t narcs = data_salbp->arcs.size();
    IloNum upper, lower;
    int number_taks_machine;
    IloAlgorithm::Status status_ret;
    bool column_added=false;
    bool timeout=false;
    for (size_t k=0; k<m; ++k) {
        //std::cout<<"load subproblems.. "<<k<<std::endl;
        pricing_problem& p =pricing_problems[k];
        number_taks_machine = data_salbp->tasks_in_machine[position];
        
        cplex_pricing.clearModel();
        cplex_pricing.extract(pricing_problems[k].model2);
        //if(k==5)cplex_pricing.exportModel("submodel.lp");
        cplex_pricing.resetTime();
        double limit_t =TIME_LIMIT-double( clock() - t_u ) / double( CLOCKS_PER_SEC );
        limit_t = limit_t>1 ? limit_t: 1;
        cplex_pricing.setParam(IloCplex::Param::TimeLimit,limit_t  );
        if(verify_timeout(timeout, t_u )){
            current_ub = UB+1;
            break;
        }
        //std::cout<<"solving princing.."<<std::endl;
        cplex_pricing.solve();
        status_ret = cplex_pricing.getStatus();
        upper = (status_ret == IloAlgorithm::Optimal)? cplex_pricing.getObjValue() : cplex_pricing.getBestObjValue();
        lower = ((status_ret == IloAlgorithm::Feasible)||
                 (status_ret == IloAlgorithm::Optimal))? cplex_pricing.getObjValue() : -1;
        current_ub += upper;
        if(verify_timeout(timeout, t_u )){
            current_ub = UB+1;
            break;
        }
        if(lower  > EPSILON){
            //std::cout<<"k "<<k+1<<" cplex_pricing: "<<cplex_pricing.getObjValue()<<" cputime: "<<cplex_pricing.getCplexTime()<<" "<<upper<<std::endl;

            translate_column( k,  position, data_salbp->qj,
                             data_salbp->arcs, data_salbp->tasks_in_machine);
            column_added=true;
        }
        position += number_taks_machine+1;
        
    }
    
    return column_added;
}

//================================================================

bool
Solver::find_feas_basis(){
    
    bool timeout;
    clock_t t_0 = clock();
    double mp_ojb=0.0;
    double current_ub;
    while(1){
        //std::cout<<"solve master.."<<std::endl;
        if(verify_timeout(timeout, t_0 )){time_to_feas=TIME_LIMIT; return false;}

        cplex.solve();
        if(cplex.getStatus() != IloAlgorithm::Optimal){
            std::cout<<"MP cplex.getStatus() == "<<cplex.getStatus()<<std::endl;
            std::abort();
        }

        mp_ojb = cplex.getObjValue();
        if(mp_ojb>=0) break;
        t_u = clock(); //avoid stopping the pricing solver.
        solve_pricing(current_ub);

        if(verify_timeout(timeout, t_0 )){time_to_feas=TIME_LIMIT; return false;}

    }
    
    time_to_feas = double( clock() - t_0 ) / double( CLOCKS_PER_SEC );
    
    IloNumArray slacks_j_(env);
    IloNumArray slacks_k_(env);
    cplex.getValues(slacks_j, slacks_j_);
    cplex.getValues(slacks_k, slacks_k_);
    
    for(size_t k =0; k<m; ++k) {
        if(slacks_k_[k]>1e-10){
            std::cout<<"ERROR::Solver::find_feas_basis():: NOT FEASIBLE"<<std::endl;
            abort();
        }
        slacks_k[k].setUB(0.0);
    }
    for(size_t j =0; j<n; ++j) {
        if(slacks_j_[j]>1e-10){
            std::cout<<"ERROR::Solver::find_feas_basis():: NOT FEASIBLE"<<std::endl;
            abort();
        }
        slacks_j[j].setUB(0.0);
    }
    return true;
}


//================================================================

bool
Solver::solve_colgen(){
    
    
    t_u = clock();
    double t = 0;
    double min_ub = UB+100;
    bool column_added = false;
    bool time_out = false;
    
    
    do{
        //cplex.exportModel("model.lp");
        if(verify_timeout(time_out, t_u ))break;
        
        //std::cout<<"solve master.."<<std::endl;
        cplex.solve();
        if(cplex.getStatus() != IloAlgorithm::Optimal){
            std::cout<<"MP cplex.getStatus() == "<<cplex.getStatus()<<std::endl;
            std::abort();
        }
        #if ADD_PRECED_CONSTRAINTS == 2
        if(add_row(data_salbp->qj, data_salbp->arcs)){

            if(verify_timeout(time_out, t_u ))break;

            cplex.setParam(IloCplex::Param::RootAlgorithm, 2);
            cplex.solve();
            if(cplex.getStatus() != IloAlgorithm::Optimal){
                std::cout<<"MP2 cplex.getStatus() == "<<cplex.getStatus()<<std::endl;
                std::abort();
            }

            if(verify_timeout(time_out, t_u ))break;
            cplex.setParam(IloCplex::Param::RootAlgorithm, 1);
        }
        #endif
        double mp_ojb = cplex.getObjValue();
        double current_ub = mp_ojb;
        //std::cout<<" rho: "<<mp_ojb;//<<std::endl;
        column_added = solve_pricing(current_ub);
        if(current_ub < min_ub){
            min_ub = current_ub;
        }
        
        //std::cout<<" current_ub: "<<current_ub<<" cputime: "<<t<<" cols "<<collection->get_size()<<std::endl;
        if(fabs(current_ub-mp_ojb)<PRES) break;
        //std::cout<<"col "<<column_added<<" disturb "<<disturbed<<" perturb "<<perturb<<" e2: "<<e2<<" ecuts "<<ecuts<<" eps "<<eps<<std::endl;
        if(verify_timeout(time_out, t_u ))break;
        
    }while((column_added) && (min_ub > LB+EPSILON));
    //std::cout<<"ecuts "<<ecuts<<" e2 "<<e2<<std::endl;
    UB = min_ub<UB ? min_ub : UB;
    return time_out;
}

//================================================================

void
Solver::solve(const std::string& instance, bool add_slacks, std::vector<int> best_id){
    if (UB<=LB+EPSILON) {
        std::ofstream outfile("fileout", std::ios::app);
        outfile<<instance<<" LB "<<LB<<" UB "
        <<UB<<" time 0.0 "<<std::endl;
        outfile.close();
        outfile.open("fileout_mip", std::ios::app);
        outfile<<instance<<" nomip"<<std::endl;
        outfile.close();
        return;
    }
    //std::cout<<"create model.."<<std::endl;
    set_parameters();
    create_model(data_salbp->arcs, data_salbp->qj, data_salbp->machines, add_slacks);
    create_subproblem(data_salbp->T, data_salbp->machines, data_salbp->tasks_in_machine);
    cplex.extract(model);
    
    //cplex.exportModel("model.lp");
    size_t num_patterns_init = collection->get_size();
    bool initfeas=true; 
    if(add_slacks)initfeas = find_feas_basis();
    bool timeout = false;
    if(initfeas)  timeout = solve_colgen();
    
    
    if(timeout || !initfeas){
        std::ofstream outfile("fileout", std::ios::app);
        outfile<<instance<<" LB "<<LB<<" UB "
        <<UB<<" time "<<TIME_LIMIT<<" num_col "<<collection->get_size()<<" time_to_feas "<<time_to_feas<<std::endl;
        //outfile<<std::endl;
        outfile.close();
    }else{
        cplex.solve();
        UB = cplex.getObjValue()<UB? UB:cplex.getObjValue() ;
        std::ofstream outfile("fileout", std::ios::app);
        outfile<<instance<<" LB "<<LB<<" UB "
        <<UB<<" time "
        <<double( clock() - t_u ) / double( CLOCKS_PER_SEC )<<" num_col "<<collection->get_size()<<" time_to_feas "<<time_to_feas<<std::endl;
        //outfile<<std::endl;
        outfile.close();
    }
    double t=double( clock() - t_u ) / double( CLOCKS_PER_SEC );
    if(UB>LB+EPSILON && t<(TIME_LIMIT-1) && initfeas){
        cplex.resetTime();
        double sol = solve_mip(t, data_salbp->qj, best_id);
        std::ofstream outfile("fileout_mip", std::ios::app);
        outfile<<instance<<" mip "<<sol<<" time "<<t<<std::endl;
        outfile.close();
    }else{
        std::ofstream outfile("fileout_mip", std::ios::app);
        outfile<<instance<<" nomip"<<std::endl;
        outfile.close();
    }
    
    //spread_new_ub(data_salbp->machines, data_salbp->tasks_in_machine);
    
    
}


//================================================================
//================================================================
//================================================================

double
Solver::solve_mip(double & time,
                  const std::vector<intpair>&qj,  std::vector<int>& best_id){
    std::list<Pattern>::const_iterator it;

    #if ADD_VALID_INEQ ==1 
    cuts.endElements();
    for(size_t k =0; k<m; ++k) {
        IloExpr constraint2(env);
        constraint2+= rho - UB;
        it = collection->patterns.begin();
        std::advance(it,collection->patterns_idx[k]);
        for (size_t l = collection->patterns_idx[k]; l<collection->patterns_idx[k+1] ; ++l){
            if(it->uncertain) constraint2 -= (fmin((*it).rho, UB) - UB)*lam[(*it).id];
            ++it;
        }
    
        cuts.add((constraint2 <= 0));
        model.add(cuts[cuts.getSize()-1]);
        constraint2.end();
    }
    #endif

    c3.endElements();
    for(size_t j =0; j<n; ++j) {
        IloExpr constraint2(env);
        constraint2 += rho - UB;
        for(size_t k = qj[j].fst; k<=qj[j].snd; ++k){
            it = collection->patterns.begin();
            std::advance(it,collection->patterns_idx[k-1]);
            for (size_t l = collection->patterns_idx[k-1]; l<collection->patterns_idx[k] ; ++l){
                if(it->uncertain && (*it).hasTask(j, collection->map)){
                    constraint2 -= (fmin((*it).rho, UB) - UB)*lam[(*it).id];
                }
                ++it;
            }
        }
        c3.add((constraint2 <= 0));
        model.add(c3[c3.getSize()-1]);
        constraint2.end();
    }

    #if ADD_PRECED_CONSTRAINTS !=1
    c4.endElements();
    int i,j;
    int cont=0;
    size_t szarcs = data_salbp->arcs.size();
    for(size_t a =0; a<szarcs;++a){
        i = data_salbp->arcs[a].fst-1;
        j = data_salbp->arcs[a].snd-1;
        //std::cout<<"arc "<<a+1<<" ("<<i+1<<","<<j+1<<")"<<" intv ["<<qj[j].fst<<","<<qj[i].snd<<"] nbcons "<<qj[i].snd-qj[j].fst+1<<std::endl;
        for(size_t k = qj[j].fst; k<=qj[i].snd; ++k){
            IloExpr constraint(env);
            //std::cout<<"station "<<k<<std::endl;
            for(size_t l = k-1; l<qj[i].snd; ++l){
                it = collection->patterns.begin();
                std::advance(it,collection->patterns_idx[l]);
                for (int p = collection->patterns_idx[l]; p<collection->patterns_idx[l+1] ; ++p){
                    if(it->hasTask(i, collection->map)){constraint += lam[it->id]; }
                    ++it;
                }
            }
            for(size_t l = k-1; l<qj[j].snd; ++l){
                it = collection->patterns.begin();
                std::advance(it,collection->patterns_idx[l]);
                for (int p = collection->patterns_idx[l]; p<collection->patterns_idx[l+1] ; ++p){
                    if(it->hasTask(j, collection->map)){constraint -= lam[it->id];  }
                    ++it;
                }
            }
            
           
            c4.add((constraint <= 0));
            model.add(c4[cont]);
            constraint.end();
            ++cont;
        }
    }
    #endif

    IloNumArray hs(env,lam.getSize());
    for(size_t j =hs.getSize(); j--; )hs[j]=0;
    while (!best_id.empty()){
        hs[best_id.back()]=1;
        best_id.pop_back();
    }
    
    model.add(IloConversion(env, lam, ILOINT));
    model.add(IloConversion(env, lam0, ILOINT));
    double limit_t =TIME_LIMIT-time;
    limit_t = limit_t>1 ? limit_t: 1;
    cplex.setParam(IloCplex::Param::TimeLimit, TIME_LIMIT-time); // Time limit in seconds
    cplex.resetTime();
    //cplex.setOut(env.out());
    //cplex.addMIPStart(lam,hs, IloCplex::MIPStartSolveFixed);
    hs.end();

    time = cplex.getCplexTime();
    cplex.solve(); 
    if(cplex.getStatus() == IloAlgorithm::Infeasible){
        std::cout<<"MIP cplex.getStatus() == IloAlgorithm::Infeasible"<<std::endl;
    }else if(cplex.getStatus() == IloAlgorithm::Feasible ||
             cplex.getStatus() == IloAlgorithm::Optimal){
        time = cplex.getCplexTime()-time;
        /* IloNumArray lam_(env);
        cplex.getValues(lam,lam_);
        std::list<Pattern>::const_iterator it;
        for(size_t k =0; k<m; ++k) {
            it = collection->patterns.begin();
            std::advance(it,collection->patterns_idx[k]);
            
            for (size_t l = collection->patterns_idx[k]; l<collection->patterns_idx[k+1] ; ++l){
                if(lam_[(*it).id]>0.7){
                    std::cout<<"k "<<k+1<<" : ";
                    const std::vector<int>& v = it->v;
                    for(size_t i =v.size();i--;){
                        std::cout<<v[i]<<" ";
                    }
                    std::cout<<std::endl;
                }
                ++it;
            }
        }
        lam_.end(); */
        return cplex.getObjValue();
    }
    return -1;
}

//================================================================

void
Solver::variable_fixing(double current_ub,
                        const IloNumArray &lam_,
                        const IloNumArray& rc){
    
    size_t numbvars = lam_.getSize();
    size_t numpatt;
    for (int k = 0; k < data_salbp->m; ++k){
        std::list<Pattern>::const_iterator it = collection->patterns.begin();
        std::advance(it,collection->patterns_idx[k]);
        numpatt = (collection->patterns_idx[k+1]-collection->patterns_idx[k]);
        if(numpatt <= 1) continue;
        for (int l = 0; l< numpatt; ++l){
            const Pattern& pattern = (*it);
            if(pattern.id < numbvars){
                if(lam_[pattern.id]<EPSILON && lam[pattern.id].getUB() > 0.9){
                    if(current_ub + rc[pattern.id] < LB-EPSILON ){
                        //std::cout<<" fix to 0: "<<pattern.rho<<(pattern.uncertain? "  uncertain ": " certain ")<<" lam: "<<lam_[pattern.id]<<std::endl;
                        lam[pattern.id].setUB(0.0);
                    }
                }else if(lam_[pattern.id]>1-EPSILON && lam[pattern.id].getLB() < 0.1){
                    if(current_ub - rc[pattern.id] < LB-EPSILON ){
                        std::cout<<" fix to 1: "<<pattern.rho<<(pattern.uncertain? "  uncertain ": " certain ")<<" lam: "<<lam_[pattern.id]<<std::endl;
                        lam[pattern.id].setLB(1.0);
                    }
                }
            }
            ++it;
        }
    }
}



//================================================================

void
Solver::spread_new_ub(const std::vector<machine> & machines,
                      const std::vector<int>&  tasks_in_machine){
    std::list<Pattern>::const_iterator it;
    size_t position=0;
    int number_taks_machine;
    for (size_t k=0; k<m; ++k) {
        
        pricing_problem& p =pricing_problems[k];
        number_taks_machine = tasks_in_machine[position];
        /*for (size_t j=0; j<number_taks_machine; ++j) {
            
            p.s1[j].setLinearCoef(p.x[j], -UB);
            
        }*/
        position += number_taks_machine+1;
        
    }
    
    for(size_t j =0; j<n; ++j) {
        IloExpr rngexpr = c3[j].getExpr();
        for (IloExpr::LinearIterator rngit(rngexpr.getLinearIterator()); rngit.ok(); ++rngit) {
            if(rngit.getCoef() > UB){
                c3[j].setLinearCoef(rngit.getVar(), -UB);
            }
        }
    }
    
    
}
