//
//  data.cpp
//  
//
//  Created by Rui Shibasaki on 15/02/21.
//

#include <stdio.h>
#include "data.hpp"
#include <queue>
#define EPS 1e-6

//================================================================
//================================================================
//================================================================


bool compare_time (double i,double j) { return (i>j); }

//================================================================

void preprocess_t(const std::vector<task>& T,
                  std::vector<double>& v1,
                  double & sum, double & sum_unc, double& max){
    sum=0; sum_unc=0;
    max=0;
    int n = T.size();
    for (int i=0; i<n ; ++i){
        const task & t = T[i];
        if(t.is_unc()){
            if(t.time >max) max =t.time;
            sum_unc+=t.time;
            v1.push_back(t.time);
        }
        sum +=t.time;
    }
    std::sort(v1.begin(), v1.end(), compare_time);
}

//================================================================

double compute_theta(const std::vector<double>& v1, int p, int modu ){
    int last = floor( (modu-1)/double(p) +EPS);
    int idx;
    double sum;
    double max=0;
    for (int k=1; k<=last; ++k){
        sum=0;
        for (int i=0; i<=k; ++i) {
            idx = k*p+1-i;
            sum += v1[idx-1];
            //std::cout<<"theta idx: "<<idx<<" "<<v1[idx-1]<<std::endl;
        }
        if(sum>max) max=sum;
    }
    return max;
}

//================================================================

double perm_sumt(const std::vector<double>& v1, int last){
    double sum=0;
    for (int i=0; i<last; ++i){
        sum+=v1[i];
    }
    return sum;
}



//================================================================

double compute_delta(const std::vector<double>& v1, int xis, int k, double sumt, double tmk){
    double diff;
    diff = sumt - perm_sumt(v1,xis) - tmk;
    return diff>0 ? diff/double(k) : 0;
}

//================================================================

void BFS(int s, int n,
         const std::vector<std::vector<int> >&  succ,
         std::vector<int>& pred_succ){
    
    std::vector<bool> visited(n, false);
    std::list<int> queue;
    
    visited[s] = true;
    queue.push_back(s);
    
    
    std::vector<int>::const_iterator jj;
    int i, j;
    while(!queue.empty()){
        i = queue.front();
        queue.pop_front();
        
        for (jj = succ[i].begin(); jj != succ[i].end(); ++jj){
            j = *jj-1;
            if (!visited[j]){
                visited[j] = true;
                queue.push_back(j);
            }
            pred_succ[s*n+(j)] = 1; // s is predecessor of j
            pred_succ[(j)*n+s] = 2; // j is successor of s
        }
    }
}

//================================================================
//================================================================
//================================================================

bool
data::compute_intervals(double LB, std::vector<intpair>& qj_temp, graph_structure& gs){
    //return compute_intervals_PatersonAlbracht( LB,  qj_temp);
    //return compute_intervals_Pirogov( LB,  qj_temp);
    return compute_intervals_1rCmax( LB,qj_temp, gs);
}

//================================================================

double
data::find_makespan(double LB, const std::vector<double> & est, 
                    std::vector<std::pair<double, int> > & tasks) const{
    double cmax=0;
    size_t njobs = tasks.size();
    if(n >= 1){
        std::sort(tasks.begin(), tasks.end());
        for(size_t j=0; j<njobs;++j){
            const task & task_j = T[tasks[j].second-1];
            double tj = task_j.time + (task_j.is_unc() ? LB : 0);
            cmax = fmax(cmax, tasks[j].first) +tj;
        }
    }
    //std::cout<<"cmax "<<cmax<<std::endl;
    return cmax;
}

//================================================================

double
data::improve_est(double t, double est, double LB, const std::vector<int>& pred) const{
    bool improved = true;
    while(improved){
        improved = false;
        int q = floor(est/cap +EPS);
        double qt = q*cap;
        //std::cout<<"improve "<<est<<" t "<<t<<" qt "<<qt<<std::endl; 
        if(est > qt + EPS){
            double min = 1e30;
            for(size_t j = pred.size(); j--; ){
                const task & task_j = T[pred[j]-1];
                double tj = task_j.time + (task_j.is_unc() ? LB : 0);
                //std::cout<<"succ "<<pred[j]<<" tj "<<tj<<" "; 
                if(min > tj) min=tj;
            }
            //std::cout<<"min "<<min<<std::endl;
            if(est + EPS < qt+min){
                est = qt+min;
                improved = true;
            }
        }
        //std::cout<<"improve2 "<<est<<" t "<<t<<" floor "<<(t+est)/cap<<" "<<(ceil((t+est)/cap -EPS) -1)<<std::endl; 
        if((ceil((t+est)/cap -EPS) -1) > q){
            //std::cout<<"improved"<<std::endl;
            est = (q+1)*cap;
            improved = true;
        }
    }
    
    return est;
}

//================================================================

void
data::add_new_arc(int from, int to, graph_structure & gs){
    gs.arcs.push_back(intpair(from,to));
    gs.pred[to-1].push_back(from);
    gs.succ[from-1].push_back(to);
    --from; --to;

    gs.pred_succ[to*n+from] = 2;
    gs.pred_succ[from*n+to] = 1;
    for (size_t j = 0; j < n; ++j){
        if(gs.pred_succ[from*n+j]==2){ // j is predecessor of 'from'
            BFS(j, n, succ, pred_succ);
        }
    }
}

//================================================================
bool
data::compute_intervals_1rCmax(double LB, std::vector<intpair>& qj_temp, graph_structure& gs){
    
    compute_topologic_order(gs);
    double t=0;
    double sum_lb, sum_ub;

    std::vector<double> est(n,0);
    std::vector<std::pair<double, int> > Pj;
    std::vector<int>::iterator it;
    for (it=gs.topologic_order.begin(); it != gs.topologic_order.end(); ++it){
        int i = *it-1;
        const task & task_1 = T[i];
        t = task_1.time + (task_1.is_unc() ? LB : 0);
        for (size_t j = 0; j < n; ++j){
            if(gs.pred_succ[i*n+j]==2){
                Pj.push_back(std::make_pair(est[j], j+1));
                //std::cout<<"pred of "<<i+1<<": "<<j+1<<std::endl;
            }
        }
        est[i] = improve_est(t, find_makespan(LB, est, Pj), LB, gs.pred[i]);
        if(est[i] + t >  m*cap) return false; 
        Pj.clear();
       //std::cout<<"est of "<<i+1<<": "<<est[i]<<std::endl;

    }
    std::vector<double> estr(n,0);
    std::vector<double> lct(n,m*cap);    
    std::vector<std::pair<double, int> > Sj;
    std::vector<int>::reverse_iterator rit;
    for (rit=gs.topologic_order.rbegin(); rit != gs.topologic_order.rend(); ++rit){
        int i = *rit-1;
        const task & task_1 = T[i];
        t = task_1.time + (task_1.is_unc() ? LB : 0);
        for (size_t j = 0; j < n; ++j){
            if(gs.pred_succ[i*n+j]==1){
                Sj.push_back(std::make_pair(estr[j], j+1));
                //std::cout<<"pred of "<<i+1<<": "<<j+1<<std::endl;
            }
        }
        estr[i] = improve_est(t, find_makespan(LB, estr, Sj),LB, gs.succ[i]);
        if(estr[i] > m*cap) return false; 
        lct[i] = m*cap - estr[i];
        Sj.clear();
        
        
    }
    bool new_arc=false;
    for (it=gs.topologic_order.begin(); it != gs.topologic_order.end() && !new_arc; ++it){
        int i = *it-1;
        const task & task_i = T[i];
        double ti = task_i.time + (task_i.is_unc() ? LB : 0);

        std::vector<int>::iterator it2;
        for (it2=gs.topologic_order.begin(); it2 != gs.topologic_order.end() && !new_arc; ++it2){
            int j = *it2-1;
            if(i==j) continue;
            const task & task_j = T[j];
            double tj = task_j.time + (task_j.is_unc() ? LB : 0);
            if((gs.pred_succ[i*n+j]==0) && (est[i]+ti > lct[j]-tj + EPS)){
                add_new_arc(j+1,i+1, gs);
                new_arc = true;
                ++gs.added_arcs;
                //std::cout<<"add arc "<<j+1<<"-"<<i+1<<std::endl;
            }
        }
    }
    if(new_arc){
        if(!compute_intervals_1rCmax( LB,qj_temp, gs)) return false;
    }else{
        for (size_t i = 0; i < n; ++i){
            int l = 1+floor(est[i]/cap +EPS);
            int u = ceil(lct[i]/cap -EPS);
            //std::cout<<"task "<<i+1<<" "<<lct[i]<<std::endl;
            if(l >  qj_temp[i].fst){
                qj_temp[i].fst = l;
            }
            if(u < qj_temp[i].snd){
                qj_temp[i].snd = u;
            }
        }
    }
    return true;

}

//================================================================

bool
data::compute_intervals_PatersonAlbracht(double LB, std::vector<intpair>& qj_temp) const{
    
    double t=0;
    double sum_lb, sum_ub;
    std::vector<intpair> qj_aux(n);
    for (size_t i = 0; i < n; ++i){
        const task & task_1 = T[i];
        sum_lb = sum_ub = task_1.time + (task_1.is_unc() ? LB : 0);
        t = task_1.time + (task_1.is_unc() ? LB : 0);
        for (size_t j = 0; j < n; ++j){
            const task & task_2 = T[j];
            if(pred_succ[i*n+j]==1) sum_ub+=task_2.time + (task_2.is_unc() ? LB : 0) ;
            else if(pred_succ[i*n+j]==2){
                sum_lb+=task_2.time + (task_2.is_unc() ? LB : 0);
            }
            //std::cout<<pred_succ[i*n+j]<<" ";
        }

        qj_aux[i].fst = ceil(sum_lb/cap -EPS);
        qj_aux[i].snd = m + 1 - ceil(sum_ub/cap -EPS);
        if(qj_aux[i].snd <qj_aux[i].fst) return false;
        if(qj_aux[i].fst >  qj_temp[i].fst){
            qj_temp[i].fst = qj_aux[i].fst;
        }
        if(qj_aux[i].snd < qj_temp[i].snd){
            qj_temp[i].snd = qj_aux[i].snd;
        }
      
    }
    
    return true;

}

//================================================================

bool
data::compute_intervals_Pirogov(double LB, std::vector<intpair>& qj_temp) const{
    
    std::vector<double> infm(n,cap*m);
    std::vector<double> supm(n,0);
    std::vector<intpair> qj_aux(n);
    
    double t=0;
    double sum_lb, sum_ub;
    double min0, max0;

    for (size_t i = 0; i < n; ++i){
        const task & task_1 = T[i];
        sum_lb = sum_ub = task_1.time + (task_1.is_unc() ? LB : 0);
        t = task_1.time + (task_1.is_unc() ? LB : 0);
        max0=0;
        for (size_t j = 0; j < n; ++j){
            const task & task_2 = T[j];
            if(pred_succ[i*n+j]==1) sum_ub+=task_2.time + (task_2.is_unc() ? LB : 0) ;
            else if(pred_succ[i*n+j]==2){
                sum_lb+=task_2.time + (task_2.is_unc() ? LB : 0);
            }
            //std::cout<<pred_succ[i*n+j]<<" ";
        }
        for (size_t jj=pred[i].size(); jj--; ) {
            size_t j = pred[i][jj]-1;
            double max1 = fmax(supm[j], cap*(ceil((supm[j] + t)/cap -EPS)-1));
            //std::cout<<i+1<<" max1 "<<max1<<" "<<cap*(ceil((supm[j]+T[i].time)/cap)-1)<<" "<<supm[j]<<std::endl;
            if(max0 < max1)max0 = max1;
        }
        
        //std::cout<<std::endl;
        //std::cout<<"sum_lb: "<<sum_lb<<" sum_ub: "<<sum_ub;
        qj_aux[i].fst = ceil(sum_lb/cap -EPS);
        qj_aux[i].snd = m + 1 - ceil(sum_ub/cap -EPS);
        supm[i] = t + fmax(cap*(qj_aux[i].fst-1), (max0));
        //std::cout<<i+1<<" max0 "<<supm[i]<<std::endl;
        //std::cout<<" ["<<qj_aux[i].fst<<","<<qj_aux[i].snd<<"]"<<std::endl;*/
    }
    
    for (size_t i = n; i--; ){
        min0=cap*m;
        t = T[i].time + (T[i].is_unc() ? LB : 0) ;
        for (size_t jj=succ[i].size(); jj--; ) {
            size_t j = succ[i][jj]-1;
            double min1 = fmin(infm[j], cap*(floor((infm[j]-t)/cap +EPS)+1));
            if(min0 > min1)min0 = min1;
        }
        infm[i] = fmin((cap*qj_aux[i].snd), (min0))-t;
    }
    
    for (size_t i = 0; i < n; ++i){
        // std::cout<<"q1 ["<<qj_temp[i].fst<<", "<<qj_temp[i].snd<<"] q2 ["<<ceil(supm[i]/cap)<<", "<<1+floor(infm[i]/cap)<<"]";
        qj_aux[i].fst = ceil(supm[i]/cap - EPS);
        qj_aux[i].snd = 1+floor(infm[i]/cap +EPS);
        if(qj_aux[i].snd <qj_aux[i].fst) return false;

        if(qj_aux[i].fst >  qj_temp[i].fst){
            qj_temp[i].fst = qj_aux[i].fst;
        }
        if(qj_aux[i].snd < qj_temp[i].snd){
            qj_temp[i].snd = qj_aux[i].snd;
        }
        //std::cout<<" then q1 ["<<qj_temp[i].fst<<", "<<qj_temp[i].snd<<"]"<<std::endl;
    }
    
    return true;

}

//================================================================

void
data::compute_unique_orig_dest(){
    
    if(machines[0].uncertain && machines[m-1].uncertain) return;
    if(m>n) return;
    bool orig, dest;
    int n_pred, n_succ;
    
    for (size_t i = 0; i < n; ++i){
        dest = orig = true;
        n_succ = n_pred = 0;
        for (size_t j = 0; j < n; ++j){
            if(pred_succ[i*n+j]==1){ dest = false; ++n_succ;}
            else if(pred_succ[i*n+j]==2){ orig = false; ++n_pred;}
        }
        if(dest && (n_pred == n-1) && !machines[m-1].uncertain){
            unique_dest = i+1;
            qj[i].fst = m;
            qj[i].snd = m;
            machines[m-1].forced.push_back(i+1);
            if(n_succ>0){
                std::cout<<"ERROR PatternGen::compute_unique_orig_dest n_succ>0 and unique dest"<<std::endl;
                abort();
            }else if(orig) std::cout<<"ATTENTION PatternGen::compute_unique_orig_dest unique dest and unique origin"<<std::endl;
        }
        if(orig && (n_succ == n-1) && !machines[0].uncertain){
            unique_orig = i+1;
            qj[i].fst = 1;
            qj[i].snd = 1;
            machines[0].forced.push_back(i+1);
            if(n_pred>0){
                std::cout<<"ERROR PatternGen::compute_unique_orig_dest n_pred>0 and unique origin"<<std::endl;
                abort();
            }else if(dest) std::cout<<"ATTENTION PatternGen::compute_unique_orig_dest unique dest and unique origin"<<std::endl;
        }
    }
    //std::cout<<"unique origin: "<<unique_orig<<" unique destin: "<<unique_dest<<std::endl;
}


//================================================================
//================================================================
//================================================================

void
data::preprocess(double LB){
    compute_unique_orig_dest();
    compute_intervals(LB, qj, this->graph_struct);
    for (size_t i = 0; i < n; ++i){
        if(qj[i].fst == qj[i].snd){
            std::vector<int>::iterator it;
            std::vector<int>& vec = machines[qj[i].snd-1].forced;
            it = find (vec.begin(), vec.end(), i+1);
            if (it == vec.end())
                vec.push_back(i+1);
        }
    }
    
}

//================================================================

void
data::reprocess(double LB){
    compute_intervals(LB, qj, this->graph_struct);
    for (size_t i = 0; i < n; ++i){
        if(qj[i].fst == qj[i].snd){
            std::vector<int>::iterator it;
            std::vector<int>& vec = machines[qj[i].snd-1].forced;
            it = find (vec.begin(), vec.end(), i+1);
            if (it == vec.end())
                vec.push_back(i+1);
        }
    }
}

//================================================================
//================================================================
//================================================================

double
data::compute_ub( double UB){
    //First compute UB1, upper bound on rho1
    
    double sumt;
    double sum_unc, sumt_unc;
    double maxt_unc;
    std::vector<double> v1;

    preprocess_t( T, v1, sumt, sum_unc, maxt_unc);
    sumt_unc = sumt - sum_unc;
    double ub1;
    double ub11 = cap - maxt_unc;
    int X = floor(sumt_unc/cap +EPS);
    
    
    if(m-m_unc <= X){
        double ub12 = cap - (sumt - cap*(m-m_unc))/double(m_unc);
        ub1 = std::min(ub11, ub12);
        //std::cout<<"UB1: "<<ub1<<" "<<ub11<<" ub12 "<<ub12<<std::endl;
        
    }else{
        double theta = compute_theta(v1, m -X -1, n_unc);
        double p1 = (sumt - cap*X)/double(m-X);
        double p2 = std::max(sum_unc/double(m-X-1),theta);
        double ub13 = cap - std::min(p1, p2);
        ub1 = std::min(ub11, ub13);
        //std::cout<<"UB1: "<<ub1<<" "<<ub11<<" ub13 "<<ub13<<" X: "<<X<<std::endl;
        
    }
    //Then compute UBinf, upper bound on rho_inf
    int gamma = ceil(n_unc/double(m) -EPS);
    int xgamma = floor(n_unc/double(gamma) +EPS);
    int xis;
    double delta;
    double maxub=0;
    std::sort(v1.begin(), v1.end());
    for (int k=1; k<=xgamma; ++k) {
        xis = std::max(gamma*k, n_unc - (gamma-1)*(m-k));
        delta = compute_delta(v1, xis, k, sumt, cap*(m-k));
        double p1 = (cap - perm_sumt(v1, xis)/double(k))/double(gamma);
        double p2 = (cap - perm_sumt(v1, gamma) - delta)/double(gamma);
        double uninfk = std::min(p1,p2);
        if(maxub < uninfk) maxub = uninfk;
    }
    double UBinf = std::min(ub1, ((xgamma >= 1) ? maxub : 1e30));
    UB = UB > UBinf ? UBinf : UB ;
    return UB;
}

//================================================================

double
data::improve_ub( double UB, double LB){
    double epsilon = 1e-1;
    double delta;
    bool nosol;
    
    while (UB-LB > epsilon) {
        delta = 0.5*(UB+LB);
        std::vector<intpair>qj_temp(n);
        for (size_t i =0; i<n; ++i){
            qj_temp[i].fst=1;
            qj_temp[i].snd=m;
        }
        nosol=false;
        graph_structure gs(this->graph_struct);
        if(compute_intervals( delta, qj_temp, gs)){
            for (size_t i =0; i<n; ++i) {
                //std::cout<<i+1<<": "<<qj_temp[i].fst<<" "<<qj_temp[i].snd<<std::endl;
                if(qj_temp[i].fst > qj_temp[i].snd){
                    nosol = true;
                    break;
                }
            }
        }else{
            nosol = true;
        }

        if (nosol) {
            UB = delta;
        }else{
            LB = delta;
        }
        //std::cout<<"delta: "<<delta<<" ub: "<<UB<<" lb: "<<LB<<std::endl;
    }
    //std::cout<<"improve_ub:: "<<UB<<std::endl;
    return UB;
}


//================================================================
//================================================================
//================================================================


void 
data::compute_topologic_order(graph_structure & gs){
    std::list<int> roots;
    std::vector<bool> ordered(n, false);
    std::vector<int> nb_pred(n, 0);
    int count=0;
    for (size_t i = 0; i < n; ++i){
        if(gs.pred[i].empty()) roots.push_back(i+1);
        nb_pred[i] = gs.pred[i].size();
    }
    while(!roots.empty()){
        int i = roots.front();
        roots.pop_front();
        ordered[i-1] = true;
        gs.topologic_order[count++] = i;
        //std::cout<<"order "<<i<<std::endl;
        for (size_t it = 0; it < gs.succ[i-1].size(); ++it){
            int j = gs.succ[i-1][it];
            if(--nb_pred[j-1]==0){ 
                roots.push_back(j);
            }
        }
    }
    if(count<n){
        std::cout<<"ERRO::Topological_order:: Cycle Detected"<<std::endl;
        abort();
    }
    /* std::cout<<"topo: ";
    for (size_t i = 0; i < topologic_order.size(); ++i){
        std::cout<<topologic_order[i]<<" ";
    }
    std::cout<<std::endl; */
}

//================================================================

void
data::read_data( const char* const argv[]){
    // Open an input file stream
    std::string instance(argv[1]);
    std::ifstream f(instance);
    
    // Read the number of task
    f >> n;
    
    // Read the task times
    T.resize(n);
    topologic_order.resize(n);
    for (size_t i = 0; i < n; ++i) {
        f >> T[i].time;
        T[i].n = i+1;
        topologic_order[i]=i+1;
    }
    
    // Read the precedence constraints
    succ.resize(n);
    pred.resize(n);
    int ii, jj; char c;
    f >> ii; f >> c; f >> jj;
    while (ii != -1) {
        succ[ii-1].push_back(jj);
        pred[jj-1].push_back(ii);
        arcs.push_back(intpair(ii,jj));
        f >> ii; f >> c; f >> jj;
    }
    
    // Read a permutation of tasks
    perm_task.resize(n);
    for (size_t i = 0; i < n; ++i)
        f >> perm_task[i];
    
    // Read the capacity and the number of bins
    f >> cap; f >> m;
    // Read a permutation of bins
    perm_bin.resize(m);
    
    for (size_t i = 0; i < m; ++i){
        f >> perm_bin[i];
        
    }
    
    // Uncertain bins
    machines.resize(m);
    m_unc = ceil(m * atof(argv[3]) -EPS);
    for (size_t i = 0; i < m_unc; ++i)
        machines[perm_bin[i]-1].uncertain = true;
    
    // Uncertain tasks
    n_unc = ceil(n * atof(argv[2]) -EPS);
    // Determine uncertain tasks
    for (size_t i = 0; i < n_unc; ++i)
        T[perm_task[i]-1].set_unc();
    
    // Close file stream
    f.close();
    pred_succ.resize(n*n,0);
    qj.resize(n);
    for (size_t i = 0; i < n; ++i){
        BFS(i, n, succ, pred_succ);
        qj[i].fst = 1;
        qj[i].snd = m;
    }
}

//================================================================

void
data::make_sets_taks_machine(){
   
    tasks_in_machine.push_back(0);
    size_t last_count_position = 0;
    for (size_t k=1; k<=m; ++k) {
        for (size_t i = 0; i < n; ++i){
            const intpair& p = qj[i];
            if(p.fst <= k && k <= p.snd){
                tasks_in_machine.push_back(i+1);
                ++tasks_in_machine[last_count_position];
                //std::cout<<"machine "<<k<<" task "<<i+1<<std::endl;
            }
        }
        tasks_in_machine.push_back(0);
        last_count_position = tasks_in_machine.size()-1;
    }
}

//================================================================
bool 
data::test_sols_validity(size_t& best_id, std::vector<std::vector<int> >& solutions){
    bool equal, solution_invalid;
    bool valid_initial_sol=false;
    for (size_t sz=solutions.size(); sz--; ) {
        if (solutions[sz].empty()){ continue;}
        
        //--check if solution exists already--
        for (size_t szz=sz; szz--; ){
            if (solutions[szz].empty()){ continue;}
            equal=true;
            for (size_t j=m*n; j--; ){
                if(solutions[sz][j]!=solutions[szz][j]){
                    equal=false;
                    break;
                }
            }
            if(equal){
                if(best_id==szz) best_id=sz;
                solutions[szz].clear();
            }
        }
        //--check if solution is valid--
        solution_invalid=false;
        std::vector<std::pair<bool,size_t> > assigned(n,std::pair<bool,size_t>(false,0));
        for (size_t k = 0; k < m; ++k){
            for (size_t j = 0; j < n; ++j){
                if (solutions[sz][k*n+j]==1) {
                    if(assigned[j].first){
                        solution_invalid=true;
                        break;
                    }
                    assigned[j].first = true;
                    assigned[j].second = k;
                    if(qj[j].fst > (k+1) || (k+1)>qj[j].snd){
                        solution_invalid=true;
                        break;
                    }
                }
            }
            if(solution_invalid) break;
        }
        if(solution_invalid){
            solutions[sz].clear();
            continue;
        }
        for (size_t j = 0; j < n; ++j){
            int kj = assigned[j].second;
            if(!assigned[j].first){
                solution_invalid=true;
                break;
            }
            for (size_t ii = pred[j].size(); ii--; ){
                size_t i = pred[j][ii]-1;
                int ki = assigned[i].second;
                if(kj < ki){
                    solution_invalid=true;
                    break;
                }
            }
            if(solution_invalid) break;
        }
        if(solution_invalid){
            solutions[sz].clear();
            continue;
        }else{
            valid_initial_sol=true;
        }
    }
    return valid_initial_sol;
}


//================================================================
//================================================================
//================================================================
