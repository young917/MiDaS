# include <random>
# include <queue>
# include "hset.hpp"

#ifndef ALGORITHMFF_HPP
#define ALGORITHMFF_HPP
using namespace std;

class Algorithm_FF{
public:
    HyperGraph *graph;
    double p,q;
    vector<int> htable;
    vector<bool> check_node;
    string eval_opt, algo_opt, outputdir;
    // vector<vector<int>> tie;

    Algorithm_FF(double p, double q, string eval_opt, string algo_opt, string outputdir, HyperGraph *graph);
    HSet* run(double target_portion);
    void burn1(int ambassador, vector<int> &burned_nodes_list, double prob);
    void burn2(int ambassador, set<int> &add, double prob, int target);
};
#endif