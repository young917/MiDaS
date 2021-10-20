#include "attribute.hpp"
#include "utils.hpp"

#ifndef GRAPH_HPP
#define GRAPH_HPP

using namespace std;

class HyperGraph {
public:
    string dataname;
    string inputpath;
    int number_of_hedges;
    int number_of_nodes;
    vector< vector<int> > node2hyperedge; 
	vector< vector<int> > hyperedge2node;
    unordered_map<string, int> nodename2index;
    unordered_map<int, string> index2nodename;
    unordered_map<int, string> edgename;
    bool exist_edgename;
    Attribute *attr;

    HyperGraph(string inputpath, string dataname);
    vector<vector<int>> get_incidence_matrix();
};
#endif