#include "hset.hpp"
#include "graph.hpp"
#include "random"
#include <unordered_map>
#include <queue>
#include <map>

#ifndef EVALUATION_HPP
#define EVALUATION_HPP

using namespace std;

class Evaluation {
public:
    int number_of_hedges;
    int number_of_nodes;
    vector< vector<int> > node2hyperedge; 
	vector< vector<int> > hyperedge2node;
    vector< unordered_map<int, int>> node2node;
    string outputdir;
    
    Evaluation(vector<int> &hypergraph_masking, HyperGraph *graph, string outputdir);
    void get_clustering_coef();
    void get_path_length();
    void count_wcc();
    // singular value
    // size of weakly connected components.
    void get_global_clustering_coef();
    void get_effective_diameter(string path);
    void get_edge_density();

private:
    int bfs(vector<bool> &check, int start_node);
};
#endif