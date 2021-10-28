# include "algorithmFF.hpp"

Algorithm_FF::Algorithm_FF(double p, double q, string eval_opt, string algo_opt, string outputdir, HyperGraph *graph){
    this->p = p;
    this->q = q;
    this->eval_opt = eval_opt;
    this->algo_opt = algo_opt;
    this->outputdir = outputdir;
    this->graph = graph;
    this->htable.resize(graph->number_of_hedges, 0);
    this->check_node.resize(graph->number_of_nodes, false);
    // this->tie.resize(graph->number_of_nodes, vector<int>(graph->number_of_nodes, 0));
    // for (int h = 0; h < graph->number_of_hedges ; h++){
    //     int hsize = graph->hyperedge2node[h].size();
    //     for (int i = 0; i < hsize; i++){
    //         for(int j = i + 1 ; j < hsize; j++){
    //             int va = graph->hyperedge2node[h][i];
    //             int vb = graph->hyperedge2node[h][j];
    //             tie[va][vb] += 1;
    //         }
    //     }
    // }
}

void Algorithm_FF::burn1(int ambassador, vector<int> &burned_nodes_list, double prob){
    random_device rd;
    mt19937 gen(rd());
    geometric_distribution<> d(1 - prob);
    default_random_engine e;
    vector<bool> burned_nodes(graph->number_of_nodes);
    vector<bool> inqueue(graph->number_of_nodes);
    vector<pair<int, int>> candidate_nodes;
    vector<int> tie(graph->number_of_nodes);

    if(algo_opt.compare("ff_c") == 0){
        queue<int> q;
        q.push(ambassador);
        inqueue[ambassador] = true;
        bool exit = false;
        while (!q.empty()){
            int v = q.front();
            q.pop();
            burned_nodes[v] = true;
            burned_nodes_list.push_back(v);
            inqueue[v] = false;

            fill(tie.begin(), tie.end(), 0);
            int vdeg = (int)graph->node2hyperedge[v].size();
            for (int hidx = 0 ; hidx < vdeg ; hidx++){
                int h = graph->node2hyperedge[v][hidx];
                int hsize = (int)graph->hyperedge2node[h].size();
                for (int nvidx = 0; nvidx < hsize ; nvidx++){
                    int nv = graph->hyperedge2node[h][nvidx];
                    if (v != nv){
                        tie[nv]++;
                    }
                }
            }
            int n = d(gen);
            candidate_nodes.clear(); // fill with unvisited neighbor nodes
            for(int nv = 0 ; nv < graph->number_of_nodes ; nv++){
                if ((burned_nodes[nv] == false) && (inqueue[nv] == false)){
                    candidate_nodes.push_back(make_pair(-tie[nv], nv));
                }
            }
            sort(candidate_nodes.begin(), candidate_nodes.end());
            for( int i = 0 ; i < min(n, (int)candidate_nodes.size()) ; i++){ // visit n ~ geometric dist. neighbors
                int next_node = candidate_nodes[i].second;
                q.push(next_node);
                inqueue[next_node] = true;
            }
        }
    }
}

void Algorithm_FF::burn2(int ambassador, set<int> &add, double prob, int target_size){
    random_device rd;
    mt19937 gen(rd());
    geometric_distribution<> d(1 - prob);
    default_random_engine e;
    vector<bool> burned_nodes(graph->number_of_nodes);
    vector<bool> inqueue(graph->number_of_nodes);
    vector<pair<int, int>> candidate_nodes;
    vector<int> tie(graph->number_of_nodes);

    if(algo_opt.compare("ff_c") == 0){
        queue<int> q;
        q.push(ambassador);
        inqueue[ambassador] = true;
        bool exit = false;
        while (!q.empty()){
            int v = q.front();
            q.pop();
            burned_nodes[v] = true;
            if (check_node[v] == false){
                check_node[v] = true;
                int deg = graph->node2hyperedge[v].size();
                for (int hidx = 0 ; hidx < deg ; hidx ++){
                    int h = graph->node2hyperedge[v][hidx];
                    htable[h]--;
                    if (htable[h] == 0){
                        add.insert(h);
                        if ((int)add.size() == target_size){
                            break;
                        }
                    }
                }
            }
            if ((int)add.size() == target_size){
                break;
            }
            inqueue[v] = false;
            fill(tie.begin(), tie.end(), 0);

            int vdeg = (int)graph->node2hyperedge[v].size();
            for (int hidx = 0 ; hidx < vdeg ; hidx++){
                int h = graph->node2hyperedge[v][hidx];
                int hsize = (int)graph->hyperedge2node[h].size();
                for (int nvidx = 0; nvidx < hsize ; nvidx++){
                    int nv = graph->hyperedge2node[h][nvidx];
                    if (v != nv){
                        tie[nv]++;
                    }
                }
            }
            int n = d(gen);
            candidate_nodes.clear(); // fill with unvisited neighbor nodes
            for(int nv = 0 ; nv < graph->number_of_nodes ; nv++){
                if ((burned_nodes[nv] == false) && (inqueue[nv] == false)){
                    candidate_nodes.push_back(make_pair(-tie[nv], nv));
                }
            }
            sort(candidate_nodes.begin(), candidate_nodes.end());
            for( int i = 0 ; i < min(n, (int)candidate_nodes.size()) ; i++){ // visit n ~ geometric dist. neighbors
                int next_node = candidate_nodes[i].second;
                q.push(next_node);
                inqueue[next_node] = true;
            }
        }
    }
}

HSet* Algorithm_FF::run(double target_portion){
    int target_size = int(floor(graph->number_of_hedges * target_portion));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, graph->number_of_nodes - 1);

    cout << "Run FF" << endl;
    auto start = chrono::steady_clock::now();

    set<int> initial_state;
    HSet *sampled = new HSet(initial_state, graph, eval_opt);
    // initialize htalbe and check_node
    for (int h = 0 ; h < graph->number_of_hedges ; h++){
        htable[h] = (int)graph->hyperedge2node[h].size();
    }
    fill(check_node.begin(), check_node.end(), false);
    
    set<int> add;
    vector<int> pool;
    while((int)add.size() < target_size){
        int ambassador = dist(gen);
        pool.clear();
        burn1(ambassador, pool, p);
        for(int i = 0 ; i < (int)pool.size() ; i++){
            int next_ambassador = pool[i];
            burn2(next_ambassador, add, q, target_size);
            if ((int)add.size() == target_size){
                break;
            }
        }
    }
    auto end = chrono::steady_clock::now();
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