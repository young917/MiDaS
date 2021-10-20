#include "algorithmES.hpp"

void AlgorithmES::step_log(string postfix, HSet *sampled){
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

void AlgorithmES::update_prop(int &hprop, int cand){
    int length = (int)algo_opt.size();
    if ( (length > 3) && (algo_opt.compare(length-3 , 3, "min") == 0) ){
        if ( (hprop == -1) || (hprop > cand) ){
            hprop = cand;
        }
    }
    else if ( (length > 3) && (algo_opt.compare(length-3 , 3, "max") == 0) ){
        if ( (hprop == -1) || (hprop < cand) ){
            hprop = cand;
        }
    }
    else if ( (length > 3) && (algo_opt.compare(length-3 , 3, "avg") == 0) ){
        hprop += cand;
    }
}

void AlgorithmES::initiate(void){
    set<int> neighbors;
    if (algo_opt.compare(0, 6, "random") == 0){
        hedges.clear();
        for(int h = 0 ; h < graph->number_of_hedges ; h++){
            hedges.push_back(h);
        }
    }
    else{
        // clean data
        for (auto &item: leaf2hedges){
            auto &vec = item.second;
            vec.clear();
        }
        leaf2hedges.clear();
        hedge_tree.clear();
        max_key = -1;
        // Fine hprop
        if (algo_opt.compare(0, 8, "line_deg") == 0){
            for (int h = 0 ; h < graph->number_of_hedges ; h++){
                int hsize = (int)graph->hyperedge2node[h].size();
                for (int vidx = 0; vidx < hsize; vidx++){
                    int v = graph->hyperedge2node[h][vidx];
                    int vdeg = (int)graph->node2hyperedge[v].size();
                    for (int hidx = 0 ; hidx < vdeg ; hidx++){
                        int nh = graph->node2hyperedge[v][hidx];
                        neighbors.insert(nh);
                    }
                }
                int hprop = (int)neighbors.size();
                neighbors.clear();
                hedge2prop[h] = hprop;
                if ((max_key == -1) || (max_key < hprop)){
                    max_key = hprop;
                }
            }
        }
        else if (algo_opt.compare(0, 4, "size") == 0){
            for (int h = 0 ; h < graph->number_of_hedges ; h++){
                int hprop = (int)graph->hyperedge2node[h].size();
                hedge2prop[h] = hprop;
                if ( (max_key == -1) || (max_key < hprop) ){
                    max_key = hprop;
                }
            }
        }
        else{
            if(algo_opt.compare(0, 4, "core") == 0){
                for (int v = 0 ; v < graph->number_of_nodes ; v++){
                    int vcore = node_core[v];
                    node2prop[v] = vcore;
                }
            }
            else if(algo_opt.compare(0, 10, "global_deg")==0){
                for(int v = 0 ; v < graph->number_of_nodes ; v++){
                    int deg = (int)graph->node2hyperedge[v].size();
                    node2prop[v] = deg;
                }
            }
            for (int h = 0 ; h < graph->number_of_hedges ; h++){
                int hprop = -1;
                if ( (algo_opt_length > 3) && (algo_opt.compare(algo_opt_length - 3, 3, "avg") == 0) ){
                    hprop = 0;
                }
                int hsize = (int)graph->hyperedge2node[h].size();
                for (int vidx = 0 ; vidx < hsize ; vidx++){
                    int v = graph->hyperedge2node[h][vidx];
                    update_prop(hprop, node2prop[v]);
                }
                if ( (algo_opt_length > 3) && (algo_opt.compare(algo_opt_length - 3, 3, "avg") == 0) ){
                    hprop = (int)round((double)hprop / hsize);
                }
                hedge2prop[h] = hprop;
                if ( (max_key == -1) || (hprop > max_key)) {
                    max_key = hprop;
                }
            }
        }
        int height = (int)ceil(log2(max_key + 1));
        int whole_size = exp2(height+1);
        leaf_start = exp2(height);    
        hedge_tree.resize(whole_size, 0);
        // store
        for(int h = 0; h < graph->number_of_hedges ; h++){
            int hprop = hedge2prop[h];
            leaf2hedges[hprop].push_back(h);
        }
        for(int k = 0 ; k <= max_key ; k++){
            hedge_tree[leaf_start + k] = pow((double)k, alpha) * (double)leaf2hedges[k].size();
            assert (hedge_tree[leaf_start + k] >= 0.0);
        }
        for(int p = leaf_start - 1 ; p > 0 ; p--){
            hedge_tree[p] = hedge_tree[2 * p] + hedge_tree[2 * p + 1];
            assert (hedge_tree[p] >= 0.0);
        }
        // string fqpath = outputdir + "frequency_init.txt";
        // ofstream pqFile(fqpath.c_str());
        // pqFile << "key,init" << endl;
        // for (int k = 0 ; k <= max_key ; k++){
        //     pqFile << to_string(k) << "," << to_string((int)leaf2hedges[k].size()) << endl;
        // }
        // pqFile.close();
    }
    return;
}

int AlgorithmES::sample_hedge(void){
    std::random_device rd;
    std::mt19937 gen(rd());

    if (algo_opt.compare(0, 6, "random") == 0){
        int numedges = (int)hedges.size();
        std::uniform_int_distribution<> dist(0, numedges-1);
        int random_index = dist(gen);
        int sampled_hedge = hedges[random_index];
        // erase hyperedge
        hedges.erase(hedges.begin() + random_index);

        return sampled_hedge;
    }
    else{
        // cou << "Hedge Sampling" << endl;
        std::uniform_real_distribution<> dist(0, 1);
        int idx = 1;
        int option_length = (int)algo_opt.size();
        while (idx < leaf_start){
            if (hedge_tree[2 * idx] == 0){
                idx = 2 * idx + 1;
                continue;
            }
            else if (hedge_tree[2 * idx + 1] == 0){
                idx = 2 * idx;
                continue;
            }
            if ( (option_length > 6) && (algo_opt.compare(option_length - 6, 6, "greedy") == 0) ){
                assert (hedge_tree[2*idx + 1] > 0);
                idx = 2 * idx + 1;
            }
            else{
                double random_double = dist(gen);
                assert ((hedge_tree[2 * idx] + hedge_tree[2 * idx + 1]) > 0.0);
                double prob = (double)hedge_tree[2 * idx] / (hedge_tree[2 * idx] + hedge_tree[2 * idx + 1]);
                if (random_double < prob) idx = 2 * idx;
                else idx = 2 * idx + 1;
                assert (hedge_tree[idx] > 0.0);
            }
        }
        int final_index = idx - leaf_start;
        // cout << "Final index is ";
        // cout << final_index << endl;


        int numedges = (int)leaf2hedges[final_index].size();
        if (numedges == 0){
            cout << "Empty leaf2hedges" << endl;
            cout << to_string(hedge_tree[idx]) << " " << final_index << endl; 
        }
        assert (numedges > 0);
        std::uniform_int_distribution<> dist2(0, numedges-1);
        int random_index = dist2(gen);
        int sampled_hedge = leaf2hedges[final_index][random_index];

        // erase
        leaf2hedges[final_index].erase(leaf2hedges[final_index].begin() + random_index);
        if ((int)leaf2hedges[final_index].size() == 0){
            hedge_tree[leaf_start + final_index] = 0.0;
        }
        else{
            hedge_tree[leaf_start + final_index] = pow((double)final_index, alpha) * (double)leaf2hedges[final_index].size();
        }
        int p = leaf_start + final_index;
        // cou << "Update" << endl;
        while (p > 1){
            int parent = (int)(floor(p/2));
            hedge_tree[parent] = hedge_tree[2 * parent] + hedge_tree[2 * parent + 1];
            p = parent;
            // cout << "Negative?" << to_string(hedge_tree[parent]) << endl;
            // cout << to_string(hedge_tree[2 * parent]) << " " << to_string(hedge_tree[2 * parent + 1]) << endl;
            assert (hedge_tree[parent] >= 0.0);
        }
        // cou << "Finally updated" << endl;
        // frequency[final_index] += 1;

        return sampled_hedge;
    }
    return -1;
}

HSet* AlgorithmES::run(double target_portion, bool output){
    int target_size = int(floor(graph->number_of_hedges * target_portion));

    cout << "Run ES" << endl;
    cout << "Algorithm Option = " << algo_opt << endl;
    cout << "Alpha = " << to_string(alpha) << endl;
    HSet *sampled;
    set<int> add_set, rm_set;

    int repeat = 1;
    // if (graph->dataname.compare(0, 4, "tags") == 0) repeat = 1;
    // else if (graph->dataname.compare(0, 7, "threads") == 0) repeat = 1;
    // else if (graph->dataname.compare(0, 6, "coauth") == 0) repeat = 1;
    
    auto start = std::chrono::steady_clock::now();
    set<int> initial_state;
    initiate(); // build tree
    sampled = new HSet(initial_state, graph, eval_opt);
    add_set.clear();
    while((int)add_set.size() < target_size){
        int selected_hyperedge = sample_hedge();
        add_set.insert(selected_hyperedge);
    }
    auto end = std::chrono::steady_clock::now();
    if ((int)add_set.size() == target_size){
        sampled->change_state(add_set, "+");
        sampled->timespent = std::chrono::duration_cast<chrono::milliseconds>(end - start);
        sampled->calculate_distribution(graph, eval_opt);
        sampled->evaluate(graph->attr); 
    }
    
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