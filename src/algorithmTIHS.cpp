# include "algorithmTIHS.hpp"
int Algorithm_TIHS::sample_hedge(vector<bool> hedge_check, int remain){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, remain - 1);
    int random_index = dist(gen);
    int sampled_hedge = -1;
    int idx = 0;
    for (int h = 0 ; h < (int)hedge_check.size() ; h++){
        if (!hedge_check[h]){
            if (idx == random_index){
                sampled_hedge = h;
                break;
            }
            idx += 1;
        }
    }
    return sampled_hedge;
}

HSet* Algorithm_TIHS::run(double target_portion){
    int target_size = int(floor(graph->number_of_hedges * target_portion));
    std::random_device rd;
    std::mt19937 gen(rd());

    auto start = chrono::steady_clock::now();
    // initialize
    set<int> initial_state;
    HSet *sampled = new HSet(initial_state, graph, eval_opt);
    vector<bool> hedge_check(graph->number_of_hedges, false);

    // start sampling
    vector<bool> check_node(graph->number_of_nodes, false);
    vector<int> check_hyperedge(graph->number_of_hedges);
    for(int h = 0 ; h < graph->number_of_hedges ; h++){
        check_hyperedge[h] = (int)graph->hyperedge2node[h].size();
    }
    set<int> pool;
    while((int)pool.size() < target_size){
        int h = sample_hedge(hedge_check, graph->number_of_hedges - (int)pool.size());
        pool.insert(h);
        hedge_check[h] = true;
        if ((int)pool.size() >= target_size){
            break;
        }
        int hsize = (int)graph->hyperedge2node[h].size();
        for(int i = 0 ; i < hsize ; i++){
            int v = graph->hyperedge2node[h][i];
            if (!check_node[v]){
                int vdeg = (int)graph->node2hyperedge[v].size();
                for (int j = 0 ; j < vdeg ; j++){
                    int h_p = graph->node2hyperedge[v][j];
                    check_hyperedge[h_p] -= 1;
                    if (check_hyperedge[h_p] == 0){
                        pool.insert(h_p);
                        hedge_check[h_p] = true;
                        if ((int)pool.size() >= target_size){
                            break;
                        }
                    }
                }
                check_node[v] = true;
            }
            if ((int)pool.size() >= target_size){
                break;
            }
        }
        if ((int)pool.size() >= target_size){
            break;
        }
    }
    auto end = chrono::steady_clock::now();
    sampled->change_state(pool, "+");
    sampled->timespent = std::chrono::duration_cast<chrono::milliseconds>(end - start);
    sampled->calculate_distribution(graph, eval_opt);
    sampled->evaluate(graph->attr);
    
    cout << "[ Result ]" << endl;
    cout << "size : " << sampled->setsize << " ; " << target_size << endl;
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