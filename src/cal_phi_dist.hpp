# include <random>
# include <ctime>
# include <sys/types.h>
# include "hset.hpp"

#ifndef CALPHIDIST_HPP
#define CALPHIDIST_HPP

using namespace std;

inline void output(vector<int> node2prop, vector<int> hyperedge2prop, HyperGraph *graph, string output_path){
    vector<pair<int,int>> sorted_nodes;
    vector<bool> hyperedge_check((int)hyperedge2prop.size());
    for (int v = 0 ; v < graph->number_of_nodes ; v++){
        sorted_nodes.push_back(make_pair(node2prop[v], v));
    }
    sort(sorted_nodes.begin(), sorted_nodes.end());
    
    ofstream outputFile(output_path.c_str());
    outputFile << "degree,avg,min,max" << endl;

    int cur_deg = -1;
    double sum_phi = 0;
    double phi_num = 0;
    double max_phi = -1;
    double min_phi = -1;
    for (int v = 0 ; v < graph->number_of_nodes ; v++){
        int deg = sorted_nodes[v].first;
        int node = sorted_nodes[v].second;
        if ((cur_deg != -1) && (cur_deg != deg)){
            // output
            outputFile << to_string(cur_deg) << "," << to_string(sum_phi / phi_num) << "," << to_string(min_phi) << "," << to_string(max_phi) << endl;
        }
        cur_deg = deg;
        for (int i = 0; i < deg ; i++){
            int h = graph->node2hyperedge[node][i];
            if (hyperedge_check[h]){
                continue;
            }
            double ln_h_phi = log(hyperedge2prop[h]);
            sum_phi += ln_h_phi;
            phi_num++;
            if ((max_phi == -1) || (max_phi < ln_h_phi)){
                max_phi = ln_h_phi;
            }
            if ((min_phi == -1) || (min_phi > ln_h_phi)){
                min_phi = ln_h_phi;
            }
            hyperedge_check[h] = true;
        }
    }
    outputFile << to_string(cur_deg) << "," << to_string(sum_phi / phi_num) << "," << to_string(min_phi) << "," << to_string(max_phi) << endl;
    outputFile.close();
}

inline void cal_phi_dist(HyperGraph *graph){
    string outputpath;

    vector<int> node2prop(graph->number_of_nodes);
    vector<int> hyperedge2prop(graph->number_of_hedges);
    
    for(int v = 0 ; v < graph->number_of_nodes ; v++){
        int deg = (int)graph->node2hyperedge[v].size();
        node2prop[v] = deg;
    }
    // min
    for (int h = 0 ; h < graph->number_of_hedges ; h++){
        int hprop = -1;
        int hsize = (int)graph->hyperedge2node[h].size();
        for (int vidx = 0 ; vidx < hsize ; vidx++){
            int v = graph->hyperedge2node[h][vidx];
            if ( (hprop > node2prop[v]) || (hprop == -1) ){
                hprop = node2prop[v];
            }
        }
        hyperedge2prop[h] = hprop;
    }
    outputpath = "./results/answer_dist/" + graph->dataname + "/phi_min_analysis.txt";
    output(node2prop, hyperedge2prop, graph, outputpath);

    // max
    for (int h = 0 ; h < graph->number_of_hedges ; h++){
        int hprop = -1;
        int hsize = (int)graph->hyperedge2node[h].size();
        for (int vidx = 0 ; vidx < hsize ; vidx++){
            int v = graph->hyperedge2node[h][vidx];
            if ( (hprop < node2prop[v]) || (hprop == -1) ){
                hprop = node2prop[v];
            }
        }
        hyperedge2prop[h] = hprop;
    }
    outputpath = "./results/answer_dist/" + graph->dataname + "/phi_max_analysis.txt";
    output(node2prop, hyperedge2prop, graph, outputpath);

    // avg
    for (int h = 0 ; h < graph->number_of_hedges ; h++){
        int hprop = 0;
        int hsize = (int)graph->hyperedge2node[h].size();
        for (int vidx = 0 ; vidx < hsize ; vidx++){
            int v = graph->hyperedge2node[h][vidx];
            hprop += node2prop[v];
        }
        hprop = (int)round((double)hprop / hsize);
        hyperedge2prop[h] = hprop;
    }
    outputpath = "./results/answer_dist/" + graph->dataname + "/phi_avg_analysis.txt";
    output(node2prop, hyperedge2prop, graph, outputpath);
}

/*
inline void cal_phi_dist(HyperGraph *graph){
    string outputpath;

    unordered_map<int, long long> phi_dist;
    vector<int> node2prop(graph->number_of_nodes);
    vector<int> hyperedge2prop(graph->number_of_hedges);
    
    for(int v = 0 ; v < graph->number_of_nodes ; v++){
        int deg = (int)graph->node2hyperedge[v].size();
        node2prop[v] = deg;
    }
    min
    for (auto &item: phi_dist){
        phi_dist[item.first] = 0;
    }
    for (int h = 0 ; h < graph->number_of_hedges ; h++){
        int hprop = -1;
        int hsize = (int)graph->hyperedge2node[h].size();
        for (int vidx = 0 ; vidx < hsize ; vidx++){
            int v = graph->hyperedge2node[h][vidx];
            if ( (hprop > node2prop[v]) || (hprop == -1) ){
                hprop = node2prop[v];
            }
        }
        hyperedge2prop[h] = hprop;
        phi_dist[hprop] = phi_dist[hprop] + (1.0 / graph->number_of_hedges);
        phi_dist[hprop] = phi_dist[hprop] + 1;
    }
    outputpath = "./results/answer_dist/" + graph->dataname + "/phi_min_analysis.txt";
    ofstream min_output(outputpath.c_str());
    min_output << "phi_min,num" << endl;
    for (auto &item: phi_dist){
        min_output << to_string(item.first) << "," << to_string(item.second) << endl;
    }
    min_output.close();

    // max
    for (auto &item: phi_dist){
        phi_dist[item.first] = 0;
    }
    for (int h = 0 ; h < graph->number_of_hedges ; h++){
        int hprop = -1;
        int hsize = (int)graph->hyperedge2node[h].size();
        for (int vidx = 0 ; vidx < hsize ; vidx++){
            int v = graph->hyperedge2node[h][vidx];
            if ( (hprop < node2prop[v]) || (hprop == -1) ){
                hprop = node2prop[v];
            }
        }
        // phi_dist[hprop] = phi_dist[hprop] + (1.0 / graph->number_of_hedges);
        phi_dist[hprop] = phi_dist[hprop] + 1;
    }
    outputpath = "./results/answer_dist/" + graph->dataname + "/phi_max_dist.txt";
    ofstream max_output(outputpath.c_str());
    max_output << "phi_max,num" << endl;
    for (auto &item: phi_dist){
        max_output << to_string(item.first) << "," << to_string(item.second) << endl;
    }
    max_output.close();

    // avg
    for (auto &item: phi_dist){
        phi_dist[item.first] = 0;
    }
    for (int h = 0 ; h < graph->number_of_hedges ; h++){
        int hprop = 0;
        int hsize = (int)graph->hyperedge2node[h].size();
        for (int vidx = 0 ; vidx < hsize ; vidx++){
            int v = graph->hyperedge2node[h][vidx];
            hprop += node2prop[v];
        }
        hprop = (int)round((double)hprop / hsize);
        // phi_dist[hprop] = phi_dist[hprop] + (1.0 / graph->number_of_hedges);
        phi_dist[hprop] = phi_dist[hprop] + 1;
    }
    outputpath = "./results/answer_dist/" + graph->dataname + "/phi_avg_dist.txt";
    ofstream avg_output(outputpath.c_str());
    avg_output << "phi_avg,num" << endl;
    for (auto &item: phi_dist){
        avg_output << to_string(item.first) << "," << to_string(item.second) << endl;
    }
    avg_output.close();
}
*/
#endif