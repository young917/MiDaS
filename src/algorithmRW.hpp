# include <random>
# include <queue>
# include "hset.hpp"

#ifndef ALGORITHMRW_HPP
#define ALGORITHMRW_HPP
using namespace std;

class Algorithm_RW{
public:
    HyperGraph *graph;
    double restart;
    string eval_opt, algo_opt, outputdir;
    int givenmaxlength;
    vector<vector<int>> node2node;
    vector<int> htable;
    vector<bool> check;

    Algorithm_RW(string eval_opt, string algo_opt, string outputdir, HyperGraph *graph, double restart, int givenmaxlength);
    HSet* run(double target_portion, bool output);
    void walk(int seed_node, int max_length, set<int> &pool, int remain);
};
#endif