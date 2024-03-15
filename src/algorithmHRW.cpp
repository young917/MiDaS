# include "algorithmHRW.hpp"

Algorithm_HRW::Algorithm_HRW(string outputdir, string algo_opt, string eval_opt, HyperGraph *graph){
    this->outputdir = outputdir;
    this->algo_opt = algo_opt;
    this->eval_opt = eval_opt;
    this->graph = graph;
    
    check.resize(graph->number_of_hedges, false);
    vcheck.resize(graph->number_of_nodes, false);
    trancheck.resize(graph->number_of_nodes, false);
}

void Algorithm_HRW::noback_run(int budget){
    std::random_device rd;
    std::mt19937 gen(rd());

    // int unit_numhedge = (int)ceil((double)graph->number_of_hedges / 10.0);
    int mu_s, e_s;
    tmp.clear();
    for (int i=0 ; i<graph->number_of_nodes ; i++){
        if (!vcheck[i]){
            tmp.push_back(i);
        }
    }
    if ((int)tmp.size() == 0){
        for (int i=0 ; i<graph->number_of_nodes ; i++){
            tmp.push_back(i);
        }
    }
    std::uniform_int_distribution<> dist(0, (int)tmp.size()-1);
    mu_s = tmp[dist(gen)];

    vector<float> weights;
    float weight_total = 0.0;
    for (int hi=0 ; hi < graph->node2hyperedge[mu_s].size() ; hi++){
        int cur_e = graph->node2hyperedge[mu_s][hi];
        int cur_size = (int)graph->hyperedge2node[cur_e].size();
        float weight = max(cur_size - 1.0, 1.0);
        weights.push_back(weight);
        weight_total += weight;
    }
    for (int hi=0 ; hi<(int)weights.size(); hi++){
        weights[hi] = weights[hi] / weight_total;
    }
    std::discrete_distribution<> d(weights.begin(), weights.end());
    int e_idx = d(gen);
    e_s = graph->node2hyperedge[mu_s][e_idx];
    if (!check[e_s]){
        num_sampled_es += 1;
        sample_order.push_back(e_s);
    }

    vcheck[mu_s] = true;
    check[e_s] = true;
    int preHe = e_s;
    int num_sampled = 0;
    while (num_sampled < budget){
        // if (num_sampled % 1000 == 0){
        //     cout << "Sampling ..." <<  to_string(num_sampled) << " / " << to_string(budget) << endl;
        // }
        // sample next m_s
        string cur_path = to_string(mu_s) + "," + to_string(e_s);
        // cout << cur_path << endl;
        int he_nei_size = (int)graph->hyperedge2node[e_s].size();

        if (he_nei_size < 2){
            std::uniform_int_distribution<> dist_tmp(0, graph->number_of_nodes - 1);
            mu_s = dist_tmp(gen);
        }
        else{
            tmp.clear();
            for (int i=0; i<he_nei_size; i++){
                int v = graph->hyperedge2node[e_s][i];
                if (v != mu_s){
                    tmp.push_back(v);
                }
            }
            std::uniform_int_distribution<> dist_tmp(0, (int)tmp.size()-1);
            int mu_s_idx = dist_tmp(gen);
            mu_s = tmp[mu_s_idx];
        }
        if (((int)graph->node2hyperedge[mu_s].size() == 1) && (check[graph->node2hyperedge[mu_s][0]])){
            std::uniform_int_distribution<> dist_tmp(0, graph->number_of_nodes - 1);
            mu_s = dist_tmp(gen);
        }
        
        // sample next e_s
        tmp.clear();
        weights.clear();
        weight_total = 0.0;
        for (int hi=0 ; hi < graph->node2hyperedge[mu_s].size() ; hi++){
            int cur_e = graph->node2hyperedge[mu_s][hi];
            if (cur_e != preHe){
                int cur_size = (int)graph->hyperedge2node[cur_e].size();
                float weight = max(cur_size - 1.0, 1.0);
                weights.push_back(weight);
                weight_total += weight;
                tmp.push_back(cur_e);
            }
        }
        for (int hi=0; hi<(int)weights.size() ; hi++){
            weights[hi] = weights[hi] / weight_total;
        }
        if ((int)weights.size() == 0){
            e_s = preHe;
        }
        else{
            std::discrete_distribution<> d_tmp(weights.begin(), weights.end());
            e_idx = d_tmp(gen);
            e_s = tmp[e_idx];
        }

        preHe = e_s;
        if (!check[e_s]){
            cout << to_string(num_sampled_es) << " / " << to_string((int)ceil((double)graph->number_of_hedges * 0.9)) << endl;
            num_sampled_es++;
            sample_order.push_back(e_s);
            // add.push_back(e_s);
            // sampled->update(add, graph, "+");
            // add.clear();
        }
        vcheck[mu_s] = true;
        check[e_s] = true;
        num_sampled += 1;
        // assert(sampled->number_of_hedges == num_sampled_es);
    }
}

void Algorithm_HRW::skip_run(int budget){
    std::random_device rd;
    std::mt19937 gen(rd());

    // int unit_numhedge = (int)ceil((double)graph->number_of_hedges / 10.0);
    int mu_s, e_s;
    tmp.clear();
    for (int i=0 ; i<graph->number_of_nodes ; i++){
        if (!vcheck[i]){
            tmp.push_back(i);
        }
    }
    if ((int)tmp.size() == 0){
        for (int i=0 ; i<graph->number_of_nodes ; i++){
            tmp.push_back(i);
        }
    }
    std::uniform_int_distribution<> dist(0, (int)tmp.size()-1);
    mu_s = tmp[dist(gen)];

    vector<float> weights;
    float weight_total = 0.0;
    for (int hi=0 ; hi < graph->node2hyperedge[mu_s].size() ; hi++){
        int cur_e = graph->node2hyperedge[mu_s][hi];
        int cur_size = (int)graph->hyperedge2node[cur_e].size();
        float weight = max(cur_size - 1.0, 1.0);
        weights.push_back(weight);
        weight_total += weight;
    }
    for (int hi=0; hi<(int)weights.size() ; hi++){
        weights[hi] = weights[hi] / weight_total;
    }
    std::discrete_distribution<> d(weights.begin(), weights.end());
    int e_idx = d(gen);
    e_s = graph->node2hyperedge[mu_s][e_idx];
    if (!check[e_s]){
        num_sampled_es += 1;
        sample_order.push_back(e_s);
    }

    int preV = mu_s;
    int preHe = e_s;
    check[e_s] = true;
    vcheck[mu_s] = true;
    
    int num_sampled = 0;
    while (num_sampled < budget){
        // if (num_sampled % 1000 == 0){
        //     cout << "Sampling ..." <<  to_string(num_sampled) << " / " << to_string(budget) << endl;
        // }
        // sample next m_s
        string cur_path = to_string(mu_s) + "," + to_string(e_s);
        int nei_size = (int)tran_path[cur_path].size();
        int he_nei_size = (int)graph->hyperedge2node[e_s].size();

        if (he_nei_size < 2){
            std::uniform_int_distribution<> dist_tmp(0, graph->number_of_nodes - 1);
            mu_s = dist_tmp(gen);
        }
        else if (nei_size < he_nei_size){
            fill(trancheck.begin(), trancheck.end(), false);
            for (int i=0; i<nei_size ; i++){
                int v = tran_path[cur_path][i];
                trancheck[v] = true;
            }
            vector<int> tmp;
            for (int i=0; i<he_nei_size; i++){
                int v = graph->hyperedge2node[e_s][i];
                if (!trancheck[v]){
                    tmp.push_back(v);
                }
            }
            std::uniform_int_distribution<> dist_tmp(0, (int)tmp.size()-1);
            int mu_s_idx = dist_tmp(gen);
            mu_s = tmp[mu_s_idx];
        }
        else{
            tran_path[cur_path].clear();
            vector<int> tmp;
            for (int i=0; i<(int)graph->hyperedge2node[e_s].size() ; i++){
                if (graph->hyperedge2node[e_s][i] != mu_s){
                    tmp.push_back(graph->hyperedge2node[e_s][i]);
                }
            }
            std::uniform_int_distribution<> dist_tmp(0, (int)tmp.size()-1);
            int mu_s_idx = dist_tmp(gen);
            mu_s = tmp[mu_s_idx];
        }
        if (((int)graph->node2hyperedge[mu_s].size() == 1) && (check[graph->node2hyperedge[mu_s][0]])){
            std::uniform_int_distribution<> dist_tmp(0, graph->number_of_nodes - 1);
            mu_s = dist_tmp(gen);
        }
        
        // sample next e_s
        weights.clear();
        weight_total = 0.0;
        for (int hi=0 ; hi < graph->node2hyperedge[mu_s].size() ; hi++){
            int cur_e = graph->node2hyperedge[mu_s][hi];
            int cur_size = (int)graph->hyperedge2node[cur_e].size();
            float weight = max(cur_size - 1.0, 1.0);
            weights.push_back(weight);
            weight_total += weight;
        }
        for (int hi=0 ; hi < (int)weights.size(); hi++){
            weights[hi] = weights[hi] / weight_total;
        }
        std::discrete_distribution<> d_tmp(weights.begin(), weights.end());
        e_idx = d_tmp(gen);
        e_s = graph->node2hyperedge[mu_s][e_idx];

        cur_path = to_string(preV) + "," + to_string(preHe);
        tran_path[cur_path].push_back(mu_s);
        cur_path = to_string(mu_s) + "," + to_string(preHe);
        tran_path[cur_path].push_back(preV);
        preV = mu_s;
        preHe = e_s;

        if (!check[e_s]){
            cout << to_string(num_sampled_es) << " / " << to_string((int)ceil((double)graph->number_of_hedges * 0.9)) << endl;
            num_sampled_es++;
            sample_order.push_back(e_s);
            // add.push_back(e_s);
            // sampled->update(add, graph, "+");
            // add.clear();
        }
        check[e_s] = true;
        vcheck[mu_s] = true;
        num_sampled += 1;
        // assert(sampled->number_of_hedges == num_sampled_es);
    }
}

HSet* Algorithm_HRW::run(double target_portion){
    int budget = graph->number_of_nodes;
    int target_size = int(floor(graph->number_of_hedges * target_portion));
    
    auto start = chrono::steady_clock::now();
    set<int> initial_state;
    HSet *sampled = new HSet(initial_state, graph, eval_opt);

    while (num_sampled_es < target_size){
        cout << "Restart" << endl;
        if (algo_opt.compare("skip") == 0){
            skip_run(budget);
        }
        else if (algo_opt.compare("noback") == 0){
            noback_run(budget);
        }
    }
    auto end = chrono::steady_clock::now();
    cout << "End" << endl;

    set<int> add;
    for (int i=0; i<target_size ; i++){
        add.insert(sample_order[i]);
    }
    sampled->change_state(add, "+");
    sampled->timespent = std::chrono::duration_cast<chrono::milliseconds>(end - start);
    sampled->calculate_distribution(graph, eval_opt);
    sampled->evaluate(graph->attr);

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