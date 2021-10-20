# include "algorithmRW.hpp"

Algorithm_RW::Algorithm_RW(string eval_opt, string algo_opt, string outputdir, HyperGraph *graph, double restart, int maxlength){
    this->eval_opt = eval_opt;
    this->algo_opt = algo_opt;
    this->outputdir = outputdir;
    this->graph = graph;
    this->restart = restart;
    this->givenmaxlength = maxlength;
    this->node2node.resize(graph->number_of_nodes);

    if (algo_opt.compare(0, 4, "rw_c") == 0){
        for(int va = 0 ; va < graph->number_of_nodes ; va++){
            for(auto h : graph->node2hyperedge[va]){
                int hsize = (int)graph->hyperedge2node[h].size();
                for(int i = 0 ;  i < hsize; i++){
                    int vb = graph->hyperedge2node[h][i];
                    if(va != vb){
                        this->node2node[va].push_back(vb);
                        this->node2node[vb].push_back(va);
                    }
                }
            }
        }
    }
}

void Algorithm_RW::walk(int seed, int max_length, set<int> &pool, int remain){
    default_random_engine e;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0, 1);
    queue<int> q;
    vector<bool> visited(graph->number_of_nodes, false);
    
    int current_remain = remain;
    q.push(seed);
    int step = 0;
    while ((step < max_length) && (q.empty() == false)){
        int current = q.front();
        q.pop();
        if(!check[current]){
            check[current] = true;
            int deg = (int)graph->node2hyperedge[current].size();
            for (int hidx = 0 ; hidx < deg ; hidx++){
                int h = graph->node2hyperedge[current][hidx];
                htable[h]--;
                if (htable[h] == 0){
                    pool.insert(h);
                    current_remain--;
                    if (current_remain == 0){
                        break;
                    }
                }
            }
        }
        if (current_remain == 0){
            break;
        }
        double random_double = dist(gen);
        if(random_double < restart){
            q.push(seed);
        }
        else if ((int)node2node[current].size() != 0){
            int deg = (int)node2node[current].size();
            std::uniform_int_distribution<> dist2(0, deg-1);
            int random_index = dist2(gen);
            int next_node = node2node[current][random_index];
            q.push(next_node);
        }
        step++;
    }
}

HSet* Algorithm_RW::run(double target_portion, bool output){
    int target_size = int(floor(graph->number_of_hedges * target_portion));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist_node(0, graph->number_of_nodes-1);

    int max_length = givenmaxlength;
    if (givenmaxlength == -1){
        max_length = 100 * graph->number_of_nodes;            
    }
    else if (givenmaxlength == 1){
        max_length = graph->number_of_nodes;
    }

    auto start = chrono::steady_clock::now();
    
    set<int> initial_state;
    HSet *sampled = new HSet(initial_state, graph, eval_opt);
    htable.clear();
    check.clear();
    htable.resize(graph->number_of_hedges);
    for (int h = 0 ; h < graph->number_of_hedges ; h++){
        htable[h] = (int)graph->hyperedge2node[h].size();
    }
    check.resize(graph->number_of_nodes);
    
    vector<int> pool;
    set<int> add;
    int remain;
    while((int)add.size() < target_size){
        // cout << to_string((int)add.size()) << " / " << target_size << endl;
        string logfname = outputdir + "log.txt";
        ofstream logFile(logfname.c_str());
        logFile << to_string((int)add.size()) << "/" << to_string(target_size) << endl;
        logFile.close();

        int seed;
        seed = dist_node(gen);
        remain = target_size - (int)add.size();
        walk(seed, max_length, add, remain);
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