#include "algorithmNS.hpp"

void AlgorithmNS::step_log(string postfix, HSet *sampled){
    string writeFile = outputdir + "log_" + postfix + ".csv";
    if (!file_exist(writeFile)){
        ofstream resultFile(writeFile.c_str());
        resultFile << "eval";
        for(auto d: sampled->evaluation){
            resultFile << "," << d.first;
        }
        resultFile << endl;
    }
    ofstream resultFile(writeFile.c_str(), fstream::app);
    resultFile << "," << to_string(sampled->evaluation_value);
    for(auto d: sampled->evaluation){
        resultFile << "," << to_string(d.second);
    }
    resultFile << endl;
}

void AlgorithmNS::initiate(void){
    if (algo_opt.compare(0, 6, "random") == 0){
        nodes.clear();
        for(int v = 0 ; v < graph->number_of_nodes ; v++){
            nodes.push_back(v);
        }
    }
    else if(algo_opt.compare(0, 10, "global_deg") == 0){
        // clean data
        for(auto &item : leaf2nodes){
            auto &vec = item.second;
            vec.clear();
        }
        leaf2nodes.clear();
        node_tree.clear();
        
        // resize map
        int max_deg = 0;
        for(int v = 0; v < graph->number_of_nodes ; v++){
            int nodedegree = (int)graph->node2hyperedge[v].size();
            max_deg = max(max_deg, nodedegree);
        }
        int height = (int)ceil(log2(max_deg + 1));
        int whole_size = exp2(height+1);
        leaf_start = exp2(height);    
        node_tree.resize(whole_size, 0);
        
        // store
        for(int v = 0; v < graph->number_of_nodes ; v++){
            int nodedegree = graph->node2hyperedge[v].size();
            leaf2nodes[nodedegree].push_back(v);
        }
        for(int d = 0 ; d <= max_deg ; d++){
            if ((d == 0) || ((int)leaf2nodes[d].size() == 0)){
                node_tree[leaf_start + d] = 0.0;
            }
            else if (alpha == 0){
                node_tree[leaf_start + d] = (double)leaf2nodes[d].size();
            }
            else{
                node_tree[leaf_start + d] = pow((double)d, alpha) * (double)leaf2nodes[d].size();
            }
            assert (node_tree[leaf_start + d] >= 0.0);
        }
        for(int p = leaf_start - 1 ; p > 0 ; p--){
            node_tree[p] = node_tree[2 * p] + node_tree[2 * p + 1];
            assert (node_tree[p] >= 0.0);
        }
    }
}

int AlgorithmNS::sample_node(void){
    std::random_device rd;
    std::mt19937 gen(rd());

    if (algo_opt.compare(0, 6, "random") == 0){
        int numnodes = (int)nodes.size();
        std::uniform_int_distribution<> dist(0, numnodes-1);
        int random_index = dist(gen);
        int sampled_node = nodes[random_index];
        nodes.erase(nodes.begin() + random_index);
        return sampled_node;
    }
    else if ((algo_opt.compare(0, 4, "core") == 0) || (algo_opt.compare(0, 10, "global_deg") == 0)){
        std::uniform_real_distribution<> dist(0, 1);
        int idx = 1;
        int option_length = (int)algo_opt.length();
        while (idx < leaf_start){
            cout << idx << endl;
            assert(node_tree[idx] > 0);
            if (node_tree[2 * idx] == 0){
                idx = 2 * idx + 1;
                continue;
            }
            else if (node_tree[2 * idx + 1] == 0){
                idx = 2 * idx;
                continue;
            }
            if ( (option_length > 6) && (algo_opt.compare(option_length - 6, 6, "greedy") == 0) ){
                assert (node_tree[2*idx + 1] > 0);
                idx = 2 * idx + 1;
            }
            else{
                double random_double = dist(gen);
                assert ((node_tree[2 * idx] + node_tree[2 * idx + 1]) > 0.0);
                double prob = (double)node_tree[2 * idx] / (node_tree[2 * idx] + node_tree[2 * idx + 1]);
                if (random_double < prob) idx = 2 * idx;
                else idx = 2 * idx + 1;
                assert (node_tree[idx] > 0.0);
            }
        }
        int final_index = idx - leaf_start;

        int numnodes = (int)leaf2nodes[final_index].size();
        std::uniform_int_distribution<> dist2(0, numnodes-1);
        int random_index = dist2(gen);
        int sampled_node = leaf2nodes[final_index][random_index];

        // erase
        leaf2nodes[final_index].erase(leaf2nodes[final_index].begin() + random_index);
        if ((final_index == 0) || ((int)leaf2nodes[final_index].size() == 0)){
            node_tree[leaf_start + final_index] = 0.0;
        }
        else if (alpha == 0){
            node_tree[leaf_start + final_index] = (double)leaf2nodes[final_index].size();
        }
        else{
            node_tree[leaf_start + final_index] = pow((double)final_index, alpha) * (double)leaf2nodes[final_index].size();
        }
        assert (node_tree[leaf_start + final_index] >= 0.0);
        int p = leaf_start + final_index;
        while (p > 1){
            int parent = (int)(floor(p/2));
            node_tree[parent] = node_tree[2 * parent] + node_tree[2 * parent + 1];
            assert (node_tree[parent] >= 0.0);
            p = parent;
        }
        return sampled_node;
    }
    return -1;
}

HSet* AlgorithmNS::run(double target_portion){
    int target_hyperedge_size = int(floor(graph->number_of_hedges * target_portion));
    random_device rd;
    mt19937 g(rd());

    cout << "Algorithm Option = " << algo_opt << endl;
    
    HSet *sampled;
    set<int> initial_state;
    vector<int> htable(graph->number_of_hedges, 0); // For induced hyperedges
    vector<int> hyperedge_pool;
    set<int> add_set, rm_set;

    auto start = chrono::steady_clock::now();
    for(int h = 0 ; h < graph->number_of_hedges; h++){
        htable[h] = (int)graph->hyperedge2node[h].size();
    }
    initiate();
    sampled = new HSet(initial_state, graph, eval_opt); // empty
    
    // Sample nodes
    while((int)add_set.size() < target_hyperedge_size){
        int sampled_node = sample_node();
        hyperedge_pool.clear();
        for(auto h : graph->node2hyperedge[sampled_node]){
            hyperedge_pool.push_back(h);
        }
        shuffle(hyperedge_pool.begin(), hyperedge_pool.end(), default_random_engine(rd()));
        for (int j = 0 ; j < (int)hyperedge_pool.size() ; j++){
            int h = hyperedge_pool[j];
            htable[h]--;
            if(htable[h] == 0){
                add_set.insert(h);
            }
            if((int)add_set.size() == target_hyperedge_size){
                break;
            }
        }
    }
    auto end = chrono::steady_clock::now();
    sampled->change_state(add_set, "+");
    sampled->timespent = std::chrono::duration_cast<chrono::milliseconds>(end - start);
    sampled->calculate_distribution(graph, eval_opt);
    sampled->evaluate(graph->attr);   
    
    cout << endl;
    cout << "[ Result ]" << endl;
    for(auto e: sampled->evaluation){
        cout << e.first << ";\t";
    }
    cout << endl;
    for(auto e: sampled->evaluation){
        cout << to_string(e.second) << ";\t";
    }
    cout << endl;

    string writeFile = outputdir + "result.txt";
    ofstream resultFile(writeFile.c_str());
    resultFile << "[ Result ]" << endl;
    for(auto e: sampled->evaluation){
        resultFile << e.first << " : " << to_string(e.second) << endl;
    }
    resultFile << "Time : " << to_string(sampled->timespent.count()) << endl;
    resultFile.close();
    return sampled;
}