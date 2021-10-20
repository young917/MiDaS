#include <cmath>
#include <cassert>
#include <chrono>
#include "graph.hpp"

#ifndef HSET_HPP
#define HSET_HPP

using namespace std;

class HSet {
public:
    vector<int> hypergraph_masking;
    int setsize;
    chrono::milliseconds timespent;
    Attribute *attr;
    double evaluation_value;
    // double degree_sign; // -1 : right shift , +1 : left shift 
    vector< int > node_degree;
    map<string, double> evaluation;
    HSet(set<int> subhypergraph, HyperGraph *graph, string eval_opt);
    HSet(HSet *hset);
    HSet(string dirpath, HyperGraph *graph, string eval_opt);
    ~HSet(){
        hypergraph_masking.clear();
        node_degree.clear();
        evaluation.clear();
        free(attr);
    }
    
    bool operator==(HSet &r) const;
    vector<int> get_hyperedgeset(void);
    void evaluate(Attribute *target_attr);
    void calculate_distribution(HyperGraph *graph, string eval_opt);
    void dynamic_update_dist(set<int> deltaset, HyperGraph *graph, string sign);
    double dynamic_update_eval(set<int> deltaset, HyperGraph *graph, string sign);
    void change_state(set<int> deltaset, string sign);
    double get_Dstat(unordered_map<int,long long> &dist1, unordered_map<int,long long> &dist2, string eval_name);
    void dynamic_degree(set<int> deltah, HyperGraph *graph, string sign);
    void dynamic_size(set<int> deltah, HyperGraph *graph, string sign);
    void dynamic_pairdegree(set<int> deltah, HyperGraph *graph, string sign);
    void dynamic_intersect(set<int> deltah, HyperGraph *graph, string sign);
    double try_degree_update(set<int> deltaset, HyperGraph *graph, string sign);
    double try_size_update(set<int> deltaset, HyperGraph *graph, string sign);
    void save_as_txt(HyperGraph *graph, string outputdir, string outputname);
};
#endif