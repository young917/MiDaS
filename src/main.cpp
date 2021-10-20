#include <sys/types.h>
// #include <boost/program_options.hpp>

#include "algorithmNS.hpp"
#include "algorithmES.hpp"
#include "algorithmFF.hpp"
#include "algorithmRW.hpp"
// #include "algorithmGreedy.hpp"
#include "algorithmMH.hpp"
#include "algorithmMH2.hpp"
#include "algorithmTIHS.hpp"
// #include "algorithmAdjust.hpp"
#include "evaluation.hpp"
#include "test.hpp"
#include "adjust_time.hpp"
#include "cal_phi_dist.hpp"

using namespace std;

void save_result (string outputdir, string portion_str, string algorithm, string dataname, string algo_opt, string exchange_opt, HSet *final, HyperGraph *graph){
    cout << "Save Result" << endl;
    // evaluate final set
    Evaluation *eval = new Evaluation(final->hypergraph_masking, graph, outputdir);
    set<int> nodes;
    for (int h = 0 ;  h < graph->number_of_hedges ; h++){
        if (final->hypergraph_masking[h] == 1){
            int hsize = (int)graph->hyperedge2node[h].size();
            for (int i = 0 ; i < hsize ; i++){
                int v = graph->hyperedge2node[h][i];
                nodes.insert(v);
            }
        }
    }
    // int num_nodes = (int)nodes.size();
    // if (num_nodes < 6000){
    //     eval->get_clustering_coef();
    //     eval->get_path_length();
    // }
    eval->get_global_clustering_coef();
    eval->get_edge_density();
    eval->count_wcc();
    /* string tmp = graph->dataname + "_" + algorithm;
    if ((int)algo_opt.size() != 0){
        tmp += "_" + algo_opt;
    }
    if ((int)exchange_opt.size() != 0){
        tmp +=  + "_" + exchange_opt;
    }
    eval->get_effective_diameter(tmp);
    */

    cout << "STORE DISTRIBUTION of data in " <<  outputdir << endl;
    final->attr->to_textfile(outputdir);

    cout << "STORE SAMPLED GRAPH" << endl;
    final->save_as_txt(graph, outputdir, "sampled_graph");
    cout << endl;

    cout << "Remove Previous Results" << endl;
    string effdiampath = outputdir + "effdiameter.txt";
    remove(effdiampath.c_str());
    string svpath = outputdir + "singular_values.txt";
    remove(svpath.c_str());
    string svfullpath = outputdir + "singular_values_full.txt";
    remove(svfullpath.c_str());
    string ovpath = outputdir + "overlapness.txt";
    remove(ovpath.c_str());
    string evpath = outputdir + "entire_evaluation.txt";
    remove(evpath.c_str());
    // string nnzpath = outputdir + "nnz.txt";
    // remove(nnzpath.c_str());    

    /* cout << "STORE RESULT IN CSV" << endl;
    string csv_fname = dataname + "_" + portion_str + ".csv";
    if( !file_exist(csv_fname) ){
        ofstream resultFile(csv_fname.c_str(),  ios::out);
        resultFile << "algorithm,algo opt,exchange opt,eval avg";
        for (auto e : final->evaluation){
            resultFile << ",eval " << e.first;
        }
        resultFile << ",timestamp" << endl;
        resultFile.close();
    }
    ofstream resultFile(csv_fname.c_str(),  ios::app);
    resultFile << algorithm << "," << algo_opt << "," << exchange_opt << "," << to_string(final->evaluation_value);
    for (auto e : final->evaluation){
        resultFile << "," << to_string(e.second);
    }
    resultFile << "," << to_string(final->timespent.count()) << endl << endl;
    resultFile.close();
    */

    return;
}

int main(int argc, char* argv[]){
    string inputpath;
    string dataname;
    // algorihtm options
    string algorithm;
    double target_portion = 0.3;
    string eval_opt = "all";
    string algo_opt = "None";
    string exchange_opt = "None";
    // ES
    double alpha = 1.0;
    // TWIDDLE
    double initial_alpha = 1.0;
    double initial_stepsize = 0.3;
    double limit_alpha = 256.0;
    double limit_stepsize = 0.1;
    double update_step_ratio = 0.1;
    double evaluation_turn = 0.05;
    int search_turn = 2; // TRandom
    bool save_csv = false;
    // pso
    double reset_rate = 0;
    // adjust
    double stepsize = 1;
    string initialization_fname = "initialization.txt";
    // RW or FF options
    int maxlength = -1;
    double restart = 0.15;
    double p = 0.51;
    double q = 0.20;
    // Greedy options
    string initsetdir = "";
    int limithour = -1;
    double exchange_portion = 0;
    // MH
    bool mh_flag = false;
    int turn = 3000;
    int inner_turn = 10;
    // Flags
    bool printing = false;
    bool testing = false;
    // Test
    int num_tries = 1000;
    int repeat = -1;

    for(int i=1; i<argc ; i++){
        string input = argv[i];
        if(input.compare("--inputpath") == 0) inputpath = argv[++i];
        else if(input.compare("--dataname") == 0) dataname = argv[++i];
        else if(input.compare("--algorithm") == 0) algorithm = argv[++i];
        else if(input.compare("--target_portion") == 0) target_portion = atof(argv[++i]);
        else if(input.compare("--eval_opt") == 0) eval_opt = argv[++i];
        else if(input.compare("--algo_opt") == 0) algo_opt = argv[++i]; // (Degree/Core)_(prob/greedy)
        // es
        else if (input.compare("--alpha") == 0) alpha = atof(argv[++i]);
        // twiddle
        else if(input.compare("--limit_alpha") == 0) limit_alpha = atof(argv[++i]);
        else if(input.compare("--limit_stepsize") == 0) limit_stepsize = atof(argv[++i]);
        else if(input.compare("--initial_stepsize") == 0) initial_stepsize = atof(argv[++i]);
        else if(input.compare("--initial_alpha") == 0) initial_alpha = atof(argv[++i]);
        else if(input.compare("--update_step_ratio") == 0) update_step_ratio = atof(argv[++i]);
        else if(input.compare("--evaluation_turn") == 0) evaluation_turn = atof(argv[++i]);
        else if(input.compare("--search_turn") == 0) search_turn = atoi(argv[++i]);
        else if(input.compare("--save_csv") == 0) save_csv = true;
        // pso
        else if (input.compare("--reset_rate") == 0) reset_rate = atof(argv[++i]);
        // adjust
        else if (input.compare("--stepsize") == 0) stepsize = atof(argv[++i]);
        else if (input.compare("--initialization_fname") == 0) initialization_fname = argv[++i];
        // mh
        else if (input.compare("--mh") == 0) mh_flag = true;
        else if(input.compare("--exchange_opt") == 0) exchange_opt = argv[++i]; // When mh flag is True, ()
        else if(input.compare("--turn") == 0) turn = atoi(argv[++i]);
        else if(input.compare("--inner_turn") == 0) inner_turn = atoi(argv[++i]);
        else if(input.compare("--initsetdir") == 0) initsetdir = argv[++i];
        // rw
        else if(input.compare("--maxlength") == 0) maxlength = atoi(argv[++i]);
        else if(input.compare("--restart") == 0) restart = atof(argv[++i]);
        // ff
        else if(input.compare("--p") == 0) p = atof(argv[++i]);
        // greedy
        else if (input.compare("--limithour") == 0) limithour = atoi(argv[++i]);
        else if (input.compare("--exchange_portion") == 0) exchange_portion = atof(argv[++i]);
        // run option
        else if(input.compare("--printing") == 0) printing = true;
        else if(input.compare("--testing") == 0) testing = true; // When developing algorithm!
        // test for dynamic computation
        else if (input.compare("--num_tries") == 0) num_tries = atoi(argv[++i]);
        else if (input.compare("--repeat") == 0) repeat = atoi(argv[++i]);
    }

    HyperGraph *graph = new HyperGraph(inputpath, dataname);
    
    // mhns & mh & mhns_pair Parameter Search!!
    string outputdir = "";
    string portion_str = "";
    bool param_flag = true;
    if (algorithm.compare("answer") == 0) param_flag = false;
    if (algorithm.compare("find_skewness") == 0) param_flag = false;
    if (algorithm.compare("prepare_coreness") == 0) param_flag = false;
    if (algorithm.compare("pattern_coreness") == 0) param_flag = false;
    if (algorithm.compare("test_dynamic") == 0) param_flag = false;
    if (algorithm.compare("adjust_time") == 0) param_flag = false;
    
    vector<string> args;
    if (testing){
        outputdir = "./results/test/";
        mkdir(outputdir.c_str(), 0776);
    }
    else if(param_flag){ //" {algorithm} / {algo_opt} != None / {dataname}_{portion} / repeat"
        args.push_back(algorithm);
        if (algo_opt.compare("None") != 0){
            args.push_back(algo_opt);
        }
        stringstream stream;
        stream << std::fixed << std::setprecision(2) << target_portion;
        portion_str = stream.str();
        args.push_back(dataname + "_" + portion_str);
        if (repeat >= 0){
            args.push_back(to_string(repeat));
        }
    }

    // BaseLine algorithm
    HSet *final = NULL;
    HSet *ret, *ret_nv;
    if(algorithm.compare("ns") == 0){
        cout << "Run NS" << endl;
        if (algo_opt.compare("global_deg") == 0){
            string alpha_str;
            stringstream stream;
            stream << std::fixed << std::setprecision(4) << alpha;
            alpha_str = stream.str();
            args[1] = args[1] + "_" + alpha_str;    
        }
        outputdir = make_directory(args);
        if (file_exist(outputdir + "sampled_graph.txt")){
            cout << "Already Exist" << endl;
            return 0;
        }
        AlgorithmNS *algo = new AlgorithmNS(alpha, eval_opt, outputdir, algo_opt, graph);
        final = algo->run(target_portion, mh_flag);
        free(algo);
        algo = NULL;
    }
    else if(algorithm.compare("es") == 0){
        cout << "Run ES" << endl;
        string alpha_str;
        stringstream stream;
        stream << std::fixed << std::setprecision(4) << alpha;
        alpha_str = stream.str();
        args[1] = args[1] + "_" + alpha_str;
        outputdir = make_directory(args);
        if (file_exist(outputdir + "sampled_graph.txt")){
            cout << "Already Exist" << endl;
            return 0;
        }
        AlgorithmES *algo = new AlgorithmES(eval_opt, outputdir, algo_opt, alpha, graph);
        final = algo->run(target_portion, mh_flag);
        free(algo);
        algo = NULL;
    }
    else if(algorithm.compare("ff") == 0){
        cout << "Run FF" << endl;
        stringstream stream;
        stream << std::fixed << std::setprecision(2) << p;
        string p_str = stream.str();

        stream.str("");
        stream << std::fixed << std::setprecision(2) << q;
        string q_str = stream.str();

        args[1] = args[1] + "_" + p_str + "_" + q_str;
        outputdir = make_directory(args);
        if (file_exist(outputdir + "sampled_graph.txt")){
            cout << "Already Exist" << endl;
            return 0;
        }
        Algorithm_FF *algo = new Algorithm_FF(p, q, eval_opt, algo_opt, outputdir, graph);
        final = algo->run(target_portion, mh_flag);
        free(algo);
        algo = NULL;
    }
    else if(algorithm.compare("rw") == 0){
        cout << "Run RW" << endl;
        args[1] = args[1] + "_" + to_string(maxlength);
        outputdir = make_directory(args);
        if (file_exist(outputdir + "sampled_graph.txt")){
            cout << "Already Exist" << endl;
            return 0;
        }
        Algorithm_RW *algo = new Algorithm_RW(eval_opt, algo_opt, outputdir, graph, restart, maxlength); 
        final = algo->run(target_portion, mh_flag);
        free(algo);
        algo = NULL;
    }
    // else if (algorithm.compare("adjust") == 0){
    //     cout << "Run Adjust" << endl;
    //     stringstream stream;
    //     stream << std::fixed << std::setprecision(1) << stepsize;
    //     string stepsize_str = stream.str();
    //     stream.str("");
    //     stream << std::fixed << std::setprecision(3) << evaluation_turn;
    //     string et_str = stream.str();
    //     args[1] = args[1] + "_" + stepsize_str + "_" + et_str;
    //     outputdir = make_directory(args);
    //     cout << outputdir << endl;
    //     if (file_exist(outputdir + "sampled_graph.txt")){
    //         cout << "Already Exist" << endl;
    //         return 0;
    //     }
    //     AlgorithmAdjust *algo = new AlgorithmAdjust(eval_opt, algo_opt, initialization_fname, stepsize, limit_alpha, evaluation_turn, outputdir, graph);
    //     final = algo->run(target_portion, mh_flag);
    //     free(algo);
    //     algo = NULL;
    // }
    else if(algorithm.compare("tihs") == 0){
        cout << "Run TIHS" << endl;
        outputdir = make_directory(args);
        if (file_exist(outputdir + "sampled_graph.txt")){
            cout << "Already Exist" << endl;
            return 0;
        }
        Algorithm_TIHS *algo = new Algorithm_TIHS(eval_opt,outputdir, graph); 
        final = algo->run(target_portion, mh_flag);
        free(algo);
        algo = NULL;
    }
    // else if (algorithm.compare(0, 6, "greedy") == 0){
    //     cout << "Run Greedy" << endl;
    //     int idx = (int)args.size() - 1;
    //     args.push_back(args[idx]);
    //     if (exchange_portion > 0){
    //         stringstream stream;
    //         stream << std::fixed << std::setprecision(2) << exchange_portion;
    //         string exchange_portion_str = stream.str();
    //         args[idx] = eval_opt + "_" + exchange_portion_str;
    //     }
    //     else{
    //         args[idx] = eval_opt;
    //     }
    //     outputdir = make_directory(args);
    //     cout << outputdir << endl;
    //     if (file_exist(outputdir + "sampled_graph.txt")){
    //         cout << "Already Exist" << endl;
    //         return 0;
    //     }
    //     AlgorithmGreedy *algo = new AlgorithmGreedy(eval_opt, algo_opt, outputdir, initsetdir, graph, exchange_portion, limithour);
    //     final = algo->run(target_portion);
    //     free(algo);
    //     algo = NULL;
    // }
    else if (algorithm.compare(0, 3, "mh2") == 0){
        if ((eval_opt.compare("degree") == 0) && (algo_opt.compare("exchange")!=0)){
            args[1] = args[1];
        }
        else{
            args[1] = args[1] + "_" + eval_opt;
        }
        outputdir = make_directory(args);
        cout << outputdir + "sampled_graph.txt" << endl;
        if (file_exist(outputdir + "sampled_graph.txt")){
            cout << "Already Exist" << endl;
            return 0;
        }
        AlgorithmMH2 *algo = new AlgorithmMH2(eval_opt, outputdir, initsetdir,algo_opt, graph);
        final = algo->run(turn, target_portion);
        free(algo);
        algo = NULL;
    }
    else if (algorithm.compare(0, 2, "mh") == 0){
        if ((eval_opt.compare("degree") == 0) && (algo_opt.compare("exchange")!=0)){
            args[1] = args[1];
        }
        else{
            args[1] = args[1] + "_" + eval_opt;
        }
        outputdir = make_directory(args);
        cout << outputdir + "sampled_graph.txt" << endl;
        if (file_exist(outputdir + "sampled_graph.txt")){
            cout << "Already Exist" << endl;
            return 0;
        }
        /*
        if (algo_opt.compare("exchange") == 0){
            initsetdir = "../results/" + initsetdir + "/";
            if( !file_exist(initsetdir + "sampled_graph.txt")){
                cout << initsetdir << "doesn't exist. So Use ES Random" << endl;
                initsetdir = "../results/es/random_1.0/" + graph->dataname + "_" + portion_str;
                if (repeat > 0){
                    initsetdir += "/" + to_string(repeat);
                }
                initsetdir += "/";
            }
            cout << "Read " << initsetdir << endl;
            initset = new HSet(initsetdir, graph);
        }
        */
        AlgorithmMH *algo = new AlgorithmMH(eval_opt, outputdir, initsetdir,algo_opt, graph);
        final = algo->run(turn, target_portion);
        free(algo);
        algo = NULL;
    }

    // After initialization => MH
    /*
    if (mh_flag){
        int target_size = int(floor(graph->number_of_hedges * target_portion));
        if (final->setsize == target_size){
            save_result(outputdir, portion_str, algorithm, dataname, algo_opt, "None", final, graph);
        }
        
        HSet *initset = new HSet(final);
        free(final);
        final = NULL;
        // string exchange_opts[2] = {"plain_local_deg", "exp_local_deg"}; 
        // determine to use "plain local deg"
        args[1] = algo_opt + "_mh";
        outputdir = make_directory(args);
        
        string algo_opt = "exchange";
        AlgorithmMH *algo = new AlgorithmMH(eval_opt, outputdir, algo_opt, graph);
        final = algo->run(turn, target_portion, initset);
        free(algo);
        algo = NULL;
    }
    */

    // Other options
    if(algorithm.compare("answer") == 0){
        cout << "Get Answer Distribution" << endl;
        string answer_outputdir = "./results/answer_dist/";
        mkdir(answer_outputdir.c_str(), 0776);
        
        answer_outputdir += graph->dataname + "/";
        mkdir(answer_outputdir.c_str(), 0776);
        
        cout << "STORE DISTRIBUTION of the answer in " <<  answer_outputdir << endl;
        graph->attr->to_textfile(answer_outputdir);

        vector<int> all(graph->number_of_hedges, 1);
        Evaluation *eval = new Evaluation(all, graph, answer_outputdir);
        // if (graph->number_of_nodes < 6000){
        //     eval->get_clustering_coef();
        //     eval->get_path_length();
        // }
        eval->count_wcc();
        eval->get_global_clustering_coef();
        string tmp = graph->dataname + "_" + algorithm;
        if ((int)algo_opt.size() != 0){
            tmp += "_" + algo_opt;
        }
        if ((int)exchange_opt.size() != 0){
            tmp +=  + "_" + exchange_opt;
        }
        // eval->get_effective_diameter(tmp);
        eval->get_edge_density();
    }
    else if (algorithm.compare("find_skewness") == 0){
        vector<int> keys;
        auto dist = graph->attr->attr_list["degree"];
        double sum_dist = 0.0;
        double cumul_dist = 0.0;
        for (auto d : dist){
            sum_dist += d.second;
            if (d.second != 0){
                keys.push_back(d.first);
            }
        }
        sort(keys.begin(), keys.end());

        double degree_median = -1;
        double degree_mean = 0;
        for (int i=0 ; i < (int)keys.size() ; i++){
            int key = keys[i];
            degree_mean += (double)(key) * (dist[key]) / sum_dist;
            cumul_dist += ((double)dist[key] / sum_dist);
            if ( (degree_median == -1) && (cumul_dist >= 0.5)){
                degree_median = key;
            }
        }
        double degree_sd = 0;
        for (int i=0 ; i < (int)keys.size() ; i++){
            int key = keys[i];
            degree_sd += pow(key - degree_mean, 2) * (double)dist[key] / sum_dist;
        }
        degree_sd = pow(degree_sd, 0.5);
        double degree_skewness = 0.0;
        for(int i=0 ; i < (int)keys.size() ; i++){
            int key = keys[i];
            degree_skewness += pow((double)(key - degree_mean) / degree_sd, 3) * (double)dist[key] / sum_dist;
        }
        cout << "[" << graph->dataname << " ]" << endl;
        cout << "Median = " <<  to_string(degree_median) << endl;
        cout << "Mean = " << to_string(degree_mean) << endl;
        cout << "SD = " << to_string(degree_sd) << endl;
        cout << "Skewness = " << to_string(degree_skewness) << endl; 
    }
    else if (algorithm.compare("prepare_coreness") == 0){
        vector<int> node_core(graph->number_of_nodes, -1); // [v] : store v's coreness
        vector<int> pos(graph->number_of_nodes, -1); // [v] : position of vertex v
        vector<int> vert(graph->number_of_nodes, -1); // [position] : vertex v
        vector<bool> check_hyperedge(graph->number_of_hedges, false);
        
        cout << "Find Node Coreness" << endl;
        int md = 0;
        for (int v = 0 ; v < graph->number_of_nodes ; v++){
            int vdeg = (int)graph->node2hyperedge[v].size();
            node_core[v] = vdeg;
            if (md < vdeg){
                md = vdeg;
            }
        }
        vector<int> bin(md + 1, 0);
        for (int v = 0 ; v < graph->number_of_nodes ; v++){
            int vdeg = (int)graph->node2hyperedge[v].size();
            bin[vdeg]++; // bin[vdeg] : store number of nodes that has vdeg
        }
        int start = 0;
        for (int d = 0 ; d <= md ; d++){
            int num = bin[d];
            bin[d] = start; // bin[d] : starting index of coreness d
            start += num;
        }

        for (int v = 0 ; v < graph->number_of_nodes ; v++){
            pos[v] = bin[node_core[v]];
            vert[pos[v]] = v;
            bin[node_core[v]]++;
        }
        for (int d = md ; d > 0 ; d--){
            bin[d] = bin[d-1]; // restore starting index
        }
        bin[0] = 0;

        // check
        int previous = -1;
        for (int i = 0 ; i < graph->number_of_nodes ; i++){
            int v = vert[i];
            assert (previous <= node_core[v]);
            if (previous < node_core[v]){ // starting index
                int vdeg = (int)graph->node2hyperedge[v].size();
                if (bin[vdeg] != i){
                    cout << "Bin " << vdeg << " = " << bin[vdeg] << "  i = " << i << endl;
                }
                assert (bin[vdeg] == i);
            }
            previous = node_core[v];
        }
        
        // core decomposition
        cout << "Number of Nodes = " << graph->number_of_nodes << endl;
        for (int i = 0 ; i < graph->number_of_nodes ; i++){
            // visit vertex according to the order of pos[v]
            int v = vert[i];
            // cout << "Set Node Core " << v;
            // cout << " = " << node_core[v] << endl;
            int vdeg = (int)graph->node2hyperedge[v].size();
            for (int hidx = 0; hidx < vdeg ; hidx++){
                int h = graph->node2hyperedge[v][hidx];
                if (!check_hyperedge[h]){
                    int hsize = (int)graph->hyperedge2node[h].size();
                    for (int nvidx = 0 ; nvidx < hsize ; nvidx++){
                        int nv = graph->hyperedge2node[h][nvidx];
                        if (node_core[nv] > node_core[v]){
                            // decrease node_core[nv] by one every time one hyperedge is stripped
                            int dnv = node_core[nv];
                            int pnv = pos[nv];
                            int pw = bin[dnv];
                            int w = vert[pw];
                            if (nv != w){
                                pos[nv] = pw;
                                pos[w] = pnv;
                                vert[pnv] = w;
                                vert[pw] = nv;
                            }
                            bin[dnv]++;
                            node_core[nv]--;
                        }
                    }
                    check_hyperedge[h] = true;
                }
            }
        }
        check_hyperedge.clear();
        bin.clear();
        pos.clear();
        vert.clear();

        // store result
        string coreness_outputdir = "./results/answer_dist/";
        mkdir(coreness_outputdir.c_str(), 0776);
        
        coreness_outputdir += graph->dataname + "/";
        mkdir(coreness_outputdir.c_str(), 0776);
        
        string writeFile = coreness_outputdir + "nodecoreness.txt";
        ofstream resultFile(writeFile.c_str());
        for (int v = 0 ; v < graph->number_of_nodes ; v++){
            resultFile << graph->index2nodename[v] << "," << node_core[v] << endl;
        }
        resultFile.close();
    }
    else if (algorithm.compare("pattern_coreness") == 0){
        vector<int> node_core(graph->number_of_nodes, 0);

        string path = "./results/answer_dist/" + graph->dataname + "/";
        mkdir(path.c_str(), 0776);
        path += "nodecoreness.txt";
        ifstream coreFile(path.c_str());
        string line;
        while (getline(coreFile, line)){
            vector<string> tokens = split(line, ',');
            int v = stoi(tokens[0]);
            int vcore = stoi(tokens[1]);
            node_core[v] = vcore;
        }
        string outputpath = "./results/answer_dist/" + graph->dataname + "/node_core_deg.txt";
        ofstream outputFile(outputpath.c_str());
        outputFile << "coreness,degree" << endl;
        for (int v = 0 ; v < graph->number_of_nodes ; v++){
            int vdeg = (int)graph->node2hyperedge[v].size();
            outputFile << node_core[v] << "," << vdeg << endl;
        }
        outputFile.close();
        cout << "Save Degree & Coreness for " << dataname << endl;
        cout << endl;
    }
    else if (algorithm.compare("test_dynamic") == 0){
        test_dynamic_change(graph, num_tries);
    }
    else if (algorithm.compare("adjust_time") == 0){
        cout << "Adjust Time" << endl;
        adjust_timespent(graph, target_portion);
    }
    else if (algorithm.compare("cal_phi_dist") == 0){
        cout << "Calculate phi dist." << endl;
        cal_phi_dist(graph);
    }
    else if (!testing){
        // save intermedidate result
        int target_size = int(floor(graph->number_of_hedges * target_portion));
        if (final->setsize == target_size){
            save_result(outputdir, portion_str, algorithm, dataname, algo_opt, "None", final, graph);
        }
    }
}