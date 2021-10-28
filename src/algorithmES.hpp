#include <random>
#include <iterator>
#include <ctime>
#include <iomanip>
#include <stdio.h>
#include <cstdlib> // for std::and() and std::srand()

# include "hset.hpp"

#ifndef ALGORITHMES_HPP
#define ALGORITHMES_HPP

using namespace std;

class AlgorithmES {
public:
    int turn;
    double alpha;
    string eval_opt;
    string outputdir;
    string algo_opt;
    int algo_opt_length;
    HyperGraph *graph;

    vector<int> hedge2prop;
    vector<int> node2prop;
    vector<int> hedges;
    vector<double> hedge_tree;
    int leaf_start;
    int max_key;
    unordered_map<int, vector<int>> leaf2hedges;

    unordered_map<int, int> frequency;

    AlgorithmES(string eval_opt, string outputdir, string algo_opt, double alpha, HyperGraph *graph){
        this->eval_opt = eval_opt;
        this->outputdir = outputdir;
        this->algo_opt = algo_opt;
        this->algo_opt_length = (int)algo_opt.size();
        this->graph = graph;
        this->alpha = alpha;
        
        this->node2prop.resize(graph->number_of_nodes);
        this->hedge2prop.resize(graph->number_of_hedges);
        
        string parampath = outputdir + "parameters.txt";
        ofstream paramFile(parampath.c_str());
        paramFile << "algo opt: " << algo_opt << endl;
        paramFile << "alpha: " << to_string(alpha) << endl;
        paramFile.close();
    }
    ~AlgorithmES(){
        for (auto &item: leaf2hedges){
            auto &vec = item.second;
            vec.clear();
        }
        leaf2hedges.clear();
        hedge_tree.clear();
        hedge2prop.clear();
        node2prop.clear();
        hedges.clear();
    }

    HSet* run(double target_portion);
    void step_log(string postfix, HSet *sampled);
    void initiate(void);
    void update_prop(int &hprop, int cand);
    int sample_hedge(void);

};
#endif