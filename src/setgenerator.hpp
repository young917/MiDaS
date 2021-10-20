#include <queue>
#include "hset.hpp"

#ifndef SETGENERATOR_HPP
#define SETGENERATOR_HPP
using namespace std;

class SetGenerator {
public:
    vector<int> node_core;
    vector<int> bin;
    vector<int> pos;
    vector<int> vert;
    vector<bool> check_hyperedge;
    SetGenerator(HyperGraph *graph);
    void write_file(string outputdir, string dataname);
};
#endif