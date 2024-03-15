# include <random>
# include <queue>
# include "hset.hpp"

#ifndef ALGORITHMHRW_HPP
#define ALGORITHMHRW_HPP
using namespace std;

class Algorithm_HRW{
public:
    HyperGraph *graph;
    string outputdir, algo_opt, eval_opt;

    unordered_map<string, vector<int>> tran_path;
    vector<int> sample_order;
    vector<int> add;
    vector<int> tmp;
    vector<bool> check;
    vector<bool> vcheck;
    vector<bool> trancheck;
    int num_sampled_es = 0;

    Algorithm_HRW(string outputdir, string algo_opt, string eval_opt, HyperGraph *graph);
    HSet* run(double accruacy);
    void skip_run(int target);
    void noback_run(int target);
};
#endif