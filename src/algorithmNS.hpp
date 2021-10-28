#include <random>
#include <iterator>
#include <ctime>
#include <iomanip>
#include <stdio.h>
#include <cstdlib> // for std::and() and std::srand()

# include "hset.hpp"

#ifndef ALGORITHMNS_HPP
#define ALGORITHMNS_HPP

using namespace std;

class AlgorithmNS {
public:
    string eval_opt;
    string outputdir;
    string algo_opt;
    HyperGraph *graph;
    double alpha;

    vector<int> nodes;
    vector<double> node_tree;
    int leaf_start;
    unordered_map<int, vector<int>> leaf2nodes;

    AlgorithmNS(double alpha, string eval_opt, string outputdir, string algo_opt, HyperGraph *graph){
        this->alpha = alpha;
        this->eval_opt = eval_opt;
        this->outputdir = outputdir;
        this->algo_opt = algo_opt;
        this->graph = graph;
    }
    HSet* run(double target_portion);

    void step_log(string postfix, HSet *sampled);
    int sample_he(HSet *sampled, int flag);
    void initiate(void);
    int sample_node(void);
};
#endif