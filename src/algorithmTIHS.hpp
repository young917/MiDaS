#include <random>
#include <iterator>
#include <ctime>
#include <iomanip>
#include <stdio.h>
#include <cstdlib> // for std::and() and std::srand()
# include "hset.hpp"

#ifndef ALGORITHMTIHS_HPP
#define ALGORITHMTIHS_HPP
using namespace std;

class Algorithm_TIHS{
public:
    HyperGraph *graph;
    string eval_opt, outputdir;

    Algorithm_TIHS(string eval_opt, string outputdir, HyperGraph *graph){
        this->eval_opt = eval_opt;
        this->outputdir = outputdir;
        this->graph = graph;
    }
    HSet* run(double target_portion);
    int sample_hedge(vector<bool> hedge_check, int remain);
};
#endif