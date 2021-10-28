# include <random>
# include "hset.hpp"

#ifndef ALGORITHMMGS_HPP
#define ALGORITHMMGS_HPP

using namespace std;

class AlgorithmMGS {
public:
    string eval_opt;
    string outputdir;
    string algo_opt;
    string sample_opt;
    string initsetdir;
    HyperGraph *graph;
    vector<int> inhedges;
    vector<int> outhedges;

    AlgorithmMGS(string eval_opt, string outputdir, string initsetdir, string algo_opt, HyperGraph *graph){
        this->eval_opt = eval_opt;
        this->outputdir = outputdir;
        this->algo_opt = algo_opt;
        this->graph = graph;
        this->initsetdir = initsetdir;
    }
    ~AlgorithmMGS(){
        inhedges.clear();
        outhedges.clear();
    }

    HSet* run(int turn, double target_portion);
    // exp or plain search
    double get_acceptance_rate(double delta, double k);
    void initiate(HSet *sampled);
    int sample_he(HSet *sampled, string opt);
    void update(set<int> &add_set, set<int> &rm_set, HSet *sampled);
    void update_partial(set<int> &set, vector<int> &hedges, string opt);

    void step_log(string postfix, bool accept, HSet *sampled, double delta, double prob);
    
};
#endif
