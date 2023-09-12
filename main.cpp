

#include "data.hpp"
#include "pattern_manager.hpp"
#include "solver.hpp"
#include "heuristic.hpp"


using namespace std;


// argv[1] : input file
// argv[2] : the ratio of uncertain tasks
// argv[3] : the ratio of uncertain bins
// argv[4] : anything (optional) read file with solutions// 

//================================================================
//================================================================
//================================================================

double getBound_from_file(std::string instance, std::string input);

double read_first_solution(const string & input_file, const data & data_salbp, double & UB, size_t & best,
                           std::vector<std::vector<int> >& firstSolution);

void rename_instance_file(std::string & instance, int n_unc, int m_unc, bool prefix_only);

void print_first_solution(double firstSolutionValue, double UB,
                          const std::string & instance, const data & data_salbp,
                          const PatternCollection& collection);
//================================================================
//================================================================
//================================================================

int main (int argc, char* argv[]) {
    
    std::string instance(argv[1]);
    data data_salbp;
    data_salbp.read_data(argv);
    //------------------------------------------
    //---------- END LECTURE ----------------
    //------------------------------------------
    
    //------------------------------------------
    // DECLARATIONS AND STARTERS
    //------------------------------------------

    PatternCollection collection(data_salbp.m, data_salbp.n);
    PatternManager patternsManager(data_salbp);
    patternsManager.setCollectionPtr(&collection);

    rename_instance_file(instance, data_salbp.n_unc, data_salbp.m_unc, false);

    
    data_salbp.preprocess(patternsManager.LB); 
    
    //---------------------------------------------------
    // FIRST FEASIBLE SOLUTION
    //---------------------------------------------------
    double firstSolutionValue=-1;
    double first_bound=data_salbp.cap;
    size_t best_sol_id=0;
    std::vector<std::vector<int> > firstSolutions(10);
    
    if(argc>4){
        //std::cout<<"!!!!!!! ATTENTION !!!!!!!: READING LB FILE: be sure the input file is correct: "<<instance+".solk"<<std::endl;
        first_bound= atof(argv[4]); //patternsManager.UB;
        //firstSolutionValue = read_first_solution("../initial_solutions/"+instance+".solk", data_salbp, first_bound, best_sol_id, firstSolutions );
        //patternsManager.set_lb(firstSolutionValue);
        patternsManager.set_ub(first_bound);
        //std::cout<<instance<<" LB "<<patternsManager.LB<<" ("<<firstSolutionValue<<") UB "<<patternsManager.UB<<std::endl;
    }else{
        patternsManager.set_ub(data_salbp.compute_ub( patternsManager.UB));
        patternsManager.set_ub(data_salbp.improve_ub(patternsManager.UB, patternsManager.LB ));

        heuristicSolution heuristic_solver(patternsManager.LB, patternsManager.UB, data_salbp);
        firstSolutionValue = heuristic_solver.solve(firstSolutions);
        best_sol_id = heuristic_solver.best_sol_id;
        if(firstSolutionValue<0){
        	std::cout<<instance<<" not solved"<<std::endl;
            std::string solfile(argv[1]);
            rename_instance_file(solfile, data_salbp.n_unc, data_salbp.m_unc, true);
            solfile = "../initial_solutions/"+solfile+".solalb1";
            firstSolutionValue = read_first_solution(solfile, data_salbp, first_bound, best_sol_id, firstSolutions );
        }else std::cout<<instance<<" solved"<<std::endl;
        patternsManager.set_lb(firstSolutionValue);
    }
    data_salbp.reprocess(patternsManager.LB);
    //std::cout<<instance<<" LB "<<patternsManager.LB<<" ("<<firstSolutionValue<<") UB "<<patternsManager.UB<<std::endl;
    
    std::vector<int> best_ids;
    bool add_slacks = false;
    if(!data_salbp.test_sols_validity(best_sol_id,  firstSolutions)){
        std::cout<<"!!!!!!! ATTENTION !!!!!!!: no VALID initial feasible solution.. abort"<<std::endl;
        add_slacks=true;
    }else{
        patternsManager.solutions2Pattern(best_sol_id, firstSolutions, best_ids);
    }
    //print_first_solution(firstSolutionValue, patternsManager.UB,  instance+".sol", data_salbp, collection);
    
    Solver solver( patternsManager.LB, patternsManager.UB, data_salbp, &patternsManager);
    data_salbp.make_sets_taks_machine();
    solver.solve(instance, add_slacks, best_ids);
    return 0;


}

//================================================================
//================================================================
//================================================================

double getBound_from_file(std::string instance, std::string input){
    std::ifstream file;
    file.open(input.c_str());
    if (!file.is_open()) {
        std::cout<<"Failure to open  datafile: "<<input;
        abort();
    }
    double lb=0;
    std::string s;
    std::istringstream ss;
    
    while(getline(file,s)){
        ss.str(s);
        ss>>s;
        ss>>lb;
        if(instance.find(s)!=std::string::npos)
            break;
        ss.clear();
    }
    file.close();
    if(lb==0){
        std::cout<<"getBound_from_file()::CAUTION: Bound=0, something may be wrong : "<<input<<std::endl;
    }
    return lb;
}

//================================================================

void rename_instance_file(std::string & instance, int n_unc, int m_unc, bool prefix_only){
    size_t cpos = instance.find("/");
    while(cpos!=std::string::npos){
        instance = instance.substr(cpos+1);
        cpos = instance.find("/");
    }
    instance = instance.substr(0,instance.size()-4);
    if(prefix_only) return;
    instance += "_"+std::to_string(n_unc);
    instance += "_"+std::to_string(m_unc);
}

//================================================================


void print_first_solution(double firstSolutionValue, double UB,
                          const std::string & instance, const data & data_salbp,
                          const PatternCollection& collection){
    
    std::ofstream outfile_solution(instance);
    outfile_solution<<firstSolutionValue<<" "<<UB<<std::endl;
    for (int k = 0; k < data_salbp.m; ++k){
        //std::cout<<"patterns for machine "<<k+1<<" "<<collection.patterns_idx[k]<<" "<<collection.patterns_idx[k+1]<<std::endl;
        outfile_solution<<k+1<<" ";
        std::list<Pattern>::const_iterator it = collection.patterns.begin();
        std::advance(it,collection.patterns_idx[k]);
        if(collection.patterns_idx[k+1]-collection.patterns_idx[k]==0){
            outfile_solution<<"0 rho: "<<UB<<" certain "<<data_salbp.cap<<std::endl;
        }
        for (int l = collection.patterns_idx[k]; l<collection.patterns_idx[k+1] ; ++l){
            const Pattern& pattern = (*it);
            double sum=0;
            outfile_solution<<pattern.v.size()<<" ";
            for(int j=pattern.v.size();j--;){
                //std::cout<<pattern.v[j]<<" ";
                outfile_solution<<pattern.v[j]<<" ";
                sum += data_salbp.T[pattern.v[j]-1].time;
            }
            //std::cout<<(pattern.uncertain? "    (uncertain ": "(certain ")<<"rho: "<<pattern.rho<<")"<<std::endl;
            //std::cout<<" rho: "<<pattern.rho<<(pattern.uncertain? "  uncertain ": " certain ")<<(pattern.uncertain ? 0 : (data_salbp.cap-sum) )<<std::endl;
            outfile_solution<<" rho: "<<pattern.rho<<(pattern.uncertain? "  uncertain ": " certain ")<<(pattern.uncertain ? 0 : (data_salbp.cap-sum) )<<std::endl;
            ++it;
        }
    }
    
    outfile_solution.close();
}

//================================================================


double
read_first_solution(const string & input_file, const data & data_salbp, double & UB, size_t & best,
                    std::vector<std::vector<int> >& firstSolutions){
    std::ifstream file;
    file.open(input_file.c_str());
    if (!file.is_open()) {
        std::cout<<"Failure to open  datafile: "<<input_file;
        abort();
    }
    
    double solution_value;
    size_t nsols;
    size_t k, ntasks, task;
    size_t m = data_salbp.m;
    size_t n = data_salbp.n;

    std::string s;
    std::istringstream ss;
    
    getline(file,s);
    ss.str(s);
    ss>>solution_value;
    ss>>UB;
    ss>>nsols;
    ss.clear();

    firstSolutions.resize(nsols);
    double best_val=-1;
    for (; nsols--; ) {
        double sol_val=1e30;
        firstSolutions[nsols].assign(n*m,0);
        while(getline(file,s)){
            ss.str(s);
            ss>>k;
            ss>>ntasks;
            --k;
            double pat_val=data_salbp.cap;
            int unc_n=0;
            for(int i =0; i<ntasks; ++i){
                ss>>task;
                --task;
                firstSolutions[nsols][k*n+task]=1;
                pat_val -= data_salbp.T[task].time;
                if(data_salbp.machines[k].uncertain || data_salbp.T[task].is_unc())
                    ++unc_n;
            }
            
            pat_val = unc_n > 0 ? pat_val/double(unc_n): UB;
            sol_val =  pat_val<sol_val ? pat_val : sol_val;
            ss.clear();
            if(k+1 == m) break;
        }
        if (sol_val > best_val) {
            best_val = sol_val;
            best = nsols;
        }
    }
    file.close();
    firstSolutions[0] = firstSolutions[best];
    firstSolutions.resize(1);
    best =0;
    return best_val;
}
