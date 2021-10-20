#include <random>
#include <iterator>
#include <ctime>
#include <iomanip>
#include <stdio.h>
#include <cstdlib> // for std::and() and std::srand()

# include "setgenerator.hpp"

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

    vector<int> node_core;
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
        
        this->node_core.resize(graph->number_of_nodes, 0);
        string path = "./results/answer_dist/" + graph->dataname + "/nodecoreness.txt";
        ifstream coreFile(path.c_str());
        string line;
        while (getline(coreFile, line)){
            vector<string> tokens = split(line, ',');
            string vname = tokens[0];
            int vindex = graph->nodename2index[vname];
            int vcore = stoi(tokens[1]);
            node_core[vindex] = vcore;
        }
    }
    HSet* run(double target_portion, bool output);

    void step_log(string postfix, HSet *sampled);
    int sample_he(HSet *sampled, int flag);
    void initiate(void);
    int sample_node(void);
};
#endif