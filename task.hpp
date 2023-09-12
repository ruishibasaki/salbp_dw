#ifndef task_hpp
#define task_hpp


class task {

public:
    
    // Parameters
    int   n;     // The number of the task
    int   time;  // The time of the task    
    bool  unc;   // "TRUE" if the task is uncertain, "FALSE" otherwise
     
	// INIT and GET methods
		  task     () {unc = false;}	
         ~task     () {};
    
    int   get_time () const {return time;}   // Return the time of the task
    void  set_unc  ()  {unc = true;}    // Set the task as uncertain
    bool  is_unc   () const {return unc;}    // Return the status of uncertainty of the task
    
    task& operator=(const task & t){
        n = t.n;
        time = t.time;
        unc = t.unc;
        return *this;
    }
};



#endif /* task_hpp */
