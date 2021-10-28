# include "algorithmMGS.hpp"

double AlgorithmMGS::get_acceptance_rate(double delta, double k){
    double prob;
    if (k < 0){
        if (delta > 0){
            prob = 1.0;
        }
        else{
            prob = 0.0;
        }
    }
    else{
        double tmp = delta * k * 1;
        if (tmp > 0) tmp = 0;
        else if (tmp < -100) tmp = -100;
        prob = exp(tmp);
    }
    
    return prob;
}

void AlgorithmMGS::step_log(string postfix, bool accept, HSet *sampled, double delta, double prob){
    string writeFile = outputdir + "log_" + postfix + ".csv";
    if (!file_exist(writeFile)){
        cout << "Write Log " << writeFile << endl;
        ofstream resultFile(writeFile.c_str());
        resultFile << "accept,size,delta,prob";
        for(auto d: sampled->evaluation){
            resultFile << "," << d.first;
        }
        resultFile << endl;
    }
    ofstream resultFile(writeFile.c_str(), fstream::app);
    if (accept){
        resultFile << "1,";
    }
    else{
        resultFile << "0,";
    }
    resultFile << to_string(sampled->setsize) << "," << to_string(delta) << ",";
    resultFile << to_string(prob) << ",";
    for(auto d: sampled->evaluation){
        resultFile << "," << to_string(d.second);
    }
    resultFile << endl;
}

void AlgorithmMGS::initiate(HSet *sampled){
    inhedges.clear();
    outhedges.clear();
    vector<int> tmp;
    
    for (int h = 0 ; h < graph->number_of_hedges ; h++){
        if (sampled->hypergraph_masking[h] == 1){
            inhedges.push_back(h);
        }
        else{
            outhedges.push_back(h);
        }
    }
}

int AlgorithmMGS::sample_he(HSet *sampled, string opt){
    int sampled_hedge = -1;
    std::random_device rd;
    std::mt19937 g(rd());

    if (opt.compare("add") == 0){
        int hsize = (int)outhedges.size();
        std::uniform_int_distribution<> dist(0, hsize-1);
        int random_index = dist(g);
        sampled_hedge = outhedges[random_index];
        // std::shuffle(outhedges.begin(), outhedges.end(), g);
        // sampled_hedge = outhedges[0];
    }
    else if (opt.compare("rm") == 0){
        int hsize = (int)inhedges.size();
        std::uniform_int_distribution<> dist(0, hsize-1);
        int random_index = dist(g);
        sampled_hedge = inhedges[random_index];
        // std::shuffle(inhedges.begin(), inhedges.end(), g);
        // sampled_hedge = inhedges[0];
    }
    return sampled_hedge;
}
void AlgorithmMGS::update_partial(set<int> &hset, vector<int> &hedges, string opt){
    if (opt.compare("add") == 0){
        for (auto h : hset){
            hedges.push_back(h);
        }
    }
    else if(opt.compare("rm") == 0){
        for (auto h : hset){
            hedges.erase(remove(hedges.begin(), hedges.end(), h), hedges.end());
        }
    }
}
void AlgorithmMGS::update(set<int> &add_set, set<int> &rm_set, HSet *sampled){
    // after change sampled set
    vector<int> temp;
    if (algo_opt.compare("exchange") == 0){
        update_partial(add_set, outhedges, "rm");
        update_partial(add_set, inhedges, "add");
        update_partial(rm_set, inhedges, "rm");
        update_partial(rm_set, outhedges, "add");
    }
    else if (algo_opt.compare("add") == 0){
        update_partial(add_set, outhedges, "rm");
    }    
    else if (algo_opt.compare("remove") == 0){
        update_partial(rm_set, inhedges, "rm");
    }
}

// Run
HSet* AlgorithmMGS :: run(int turn, double target_portion){
    int target_hyperedge_size = int(floor(graph->number_of_hedges * target_portion));
    int log_step = int(floor(graph->number_of_hedges * 0.1));
    
    HSet *best = NULL;
    HSet *sampled;
    set<int> add_set, rm_set;
    int best_k, best_turn;

    cout << "\n[ MGS " + algo_opt + " ] " << endl; 

    int max_turn = turn;
    if (algo_opt.compare("add") == 0){
        max_turn = target_hyperedge_size * 1000;
    }
    else if (algo_opt.compare("remove") == 0){
        max_turn = (graph->number_of_hedges - target_hyperedge_size) * 100;
    }
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0, 1);

    double best_eval_value = 1000.0;
    int klist[4] = {10000, 1000, 100, 10};
    for (int idx = 0; idx < 4; idx ++){
        int k = klist[idx];
        cout << "[ " << k  << " ]" << endl;

        string postfix = to_string(k);
        string log_output = outputdir + "log_" + postfix + ".csv";
        int success_rm = remove(log_output.c_str());
        if (success_rm != 0){
            cout << "Fail removing the log file" << endl;
        }
        
        auto start = chrono::steady_clock::now();
        if ((int)initsetdir.size() > 0){
            sampled = new HSet(initsetdir, graph, eval_opt);
        }
        else if (algo_opt.compare("exchange") == 0){
            vector<int> initial_state_vec;
            for (int h = 0  ; h < graph->number_of_hedges ; h++){
                initial_state_vec.push_back(h);
            }
            std::shuffle(initial_state_vec.begin(), initial_state_vec.end(), gen);
            set<int> initial_state;
            for (int hi = 0 ; hi < target_hyperedge_size ; hi++){
                initial_state.insert(initial_state_vec[hi]);
            }
            sampled = new HSet(initial_state, graph, eval_opt);
        }
        else if (algo_opt.compare("add") == 0){
            set<int> initial_state;
            sampled = new HSet(initial_state, graph, eval_opt);
        }
        else if (algo_opt.compare("remove") == 0){
            set<int> initial_state;
            for (int h = 0 ; h < graph->number_of_hedges ; h++){
                initial_state.insert(h);
            }
            sampled = new HSet(initial_state, graph, eval_opt);
        }
        sampled->evaluate(graph->attr);
        if (algo_opt.compare("exchange") == 0){
            if ((best == NULL) || (best_eval_value > sampled->evaluation_value)){
                best = new HSet(sampled);
                best_eval_value = sampled->evaluation_value;
            }
        }
        initiate(sampled); // separate in_hes and out_hes

        // start exchange
        int t = 0;
        int add_h, rm_h;

        while(t < max_turn){
            if (t % 100 == 0){
                cout << t << "/" << max_turn << "      ";
            }
            if ( (algo_opt.compare("exchange") != 0) && (sampled->setsize == target_hyperedge_size)){
                break;
            }
            double previous_eval = sampled->evaluation_value;
            add_set.clear();
            rm_set.clear();
            if (algo_opt.compare("add") != 0){
                rm_h = sample_he(sampled, "rm");
                rm_set.insert(rm_h);
                sampled->dynamic_update_eval(rm_set, graph, "-");  
                // changed += sampled->try_degree_update(rm_set, graph, "-", eval_opt);  
            }
            if (algo_opt.compare("remove") != 0){
                add_h = sample_he(sampled, "add");
                add_set.insert(add_h);
                sampled->dynamic_update_eval(add_set, graph, "+"); 
                // changed += sampled->try_degree_update(add_set, graph, "+", eval_opt);   
            }
            double current_eval = sampled->evaluation_value;
            double changed = previous_eval - current_eval;
            
            double prob = get_acceptance_rate(changed, k);
            double random_double = dist(gen);
            if (random_double < prob){
                step_log(postfix, true, sampled, changed, prob);
                update(add_set, rm_set, sampled);
                cout << "Update : " << to_string(changed) << "\tsetsize=" << sampled->setsize << endl;
                for(auto e: sampled->evaluation){
                    cout << e.first << " : " << to_string(e.second) << "\t";
                }
                cout << endl;
                if ( (sampled->setsize == target_hyperedge_size) && (sampled->evaluation_value < best_eval_value)){
                    best_eval_value = sampled->evaluation_value;
                    if (best != NULL){
                        free(best);
                        best = NULL;
                    }
                    best = new HSet(sampled);
                    auto end = chrono::steady_clock::now();
                    best->timespent = chrono::duration_cast<chrono::milliseconds>(end - start);
                    best_k = k;
                    best_turn = t;
                    best->save_as_txt(graph, outputdir, "best_graph_" + to_string(best_k));
                }
            }
            else{ // restore
                step_log(postfix, false, sampled, changed, prob);
                if (algo_opt.compare("add") != 0){
                    sampled->dynamic_update_eval(rm_set, graph, "+");    
                }
                if (algo_opt.compare("remove") != 0){
                    sampled->dynamic_update_eval(add_set, graph, "-");    
                }
            }
            t++;
        }
        free(sampled);
        sampled = NULL;
    }
    if (eval_opt.compare("degree") == 0){
        best->calculate_distribution(graph, "avg");
        best->evaluate(graph->attr);
    }
    cout << endl;
    cout << "[ MGS Result ]" << endl;
    cout << "Best k : " << best_k << endl;
    cout << "Best turn : " << best_turn << endl;
    for(auto e: best->evaluation){
        cout << e.first << ";\t";
    }
    cout << endl;
    for(auto e: best->evaluation){
        cout << to_string(e.second) << ";\t";
    }
    cout << endl;

    string writeFile = outputdir + "result.txt";
    ofstream resultFile(writeFile.c_str());
    resultFile << "[ Result ]" << endl;
    resultFile << "Best k : " << best_k << endl;
    resultFile << "Best turn : " << best_turn << endl;
    for(auto e: best->evaluation){
        resultFile << e.first << " : " << to_string(e.second) << endl;
    }
    resultFile << "Time : " << to_string(best->timespent.count()) << endl;
    resultFile.close();

    return best;
}