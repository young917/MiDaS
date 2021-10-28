#include <sys/types.h>
// #include <boost/program_options.hpp>

#include "algorithmNS.hpp"
#include "algorithmES.hpp"
#include "algorithmFF.hpp"
#include "algorithmRW.hpp"
#include "algorithmMGS.hpp"
#include "algorithmTIHS.hpp"
#include "evaluation.hpp"
#include "test.hpp"
#include "adjust_time.hpp"
#include "cal_phi_dist.hpp"

using namespace std;

void save_result (string outputdir, string portion_str, string algorithm, string dataname, string algo_opt, HSet *final, HyperGraph *graph){
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
    eval->get_global_clustering_coef();
    eval->get_edge_density();
    eval->count_wcc();
    
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
    
    return;
}

int main(int argc, char* argv[]){
    string inputpath;
    string dataname;
    // algorihtm options
    string algorithm;
    double target_portion = 0.3;
    string eval_opt = "avg";
    string algo_opt = "None";
    string initsetdir = "";
    // ES
    double alpha = 1.0;
    // RW or FF options
    int maxlength = -1;
    double restart = 0.15;
    double p = 0.51;
    double q = 0.20;
    // MGS
    int turn = 3000;
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
        else if(input.compare("--algo_opt") == 0) algo_opt = argv[++i]; 
        // es
        else if (input.compare("--alpha") == 0) alpha = atof(argv[++i]);
        // mgs
        else if(input.compare("--turn") == 0) turn = atoi(argv[++i]);
        else if(input.compare("--initsetdir") == 0) initsetdir = argv[++i];
        // rw
        else if(input.compare("--maxlength") == 0) maxlength = atoi(argv[++i]);
        else if(input.compare("--restart") == 0) restart = atof(argv[++i]);
        // ff
        else if(input.compare("--p") == 0) p = atof(argv[++i]);
        else if(input.compare("--q") == 0) q = atof(argv[++i]);
        // run option
        else if(input.compare("--printing") == 0) printing = true;
        else if(input.compare("--testing") == 0) testing = true; // When developing algorithm!
        // test for dynamic computation
        else if (input.compare("--num_tries") == 0) num_tries = atoi(argv[++i]);
        else if (input.compare("--repeat") == 0) repeat = atoi(argv[++i]);
    }

    HyperGraph *graph = new HyperGraph(inputpath, dataname);
    
    string outputdir = "";
    string portion_str = "";
    
    vector<string> args;

    if (testing){
        outputdir = "./results/test/";
        mkdir(outputdir.c_str(), 0776);
    }
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
        final = algo->run(target_portion);
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
        final = algo->run(target_portion);
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
        final = algo->run(target_portion);
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
        final = algo->run(target_portion);
        free(algo);
        algo = NULL;
    }
    else if(algorithm.compare("tihs") == 0){
        cout << "Run TIHS" << endl;
        outputdir = make_directory(args);
        if (file_exist(outputdir + "sampled_graph.txt")){
            cout << "Already Exist" << endl;
            return 0;
        }
        Algorithm_TIHS *algo = new Algorithm_TIHS(eval_opt,outputdir, graph); 
        final = algo->run(target_portion);
        free(algo);
        algo = NULL;
    }
    else if (algorithm.compare(0, 3, "mgs") == 0){
        args[1] = args[1] + "_" + eval_opt;
        outputdir = make_directory(args);
        cout << outputdir + "sampled_graph.txt" << endl;
        if (file_exist(outputdir + "sampled_graph.txt")){
            cout << "Already Exist" << endl;
            return 0;
        }
        AlgorithmMGS *algo = new AlgorithmMGS(eval_opt, outputdir, initsetdir, algo_opt, graph);
        final = algo->run(turn, target_portion);
        free(algo);
        algo = NULL;
    }

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
            save_result(outputdir, portion_str, algorithm, dataname, algo_opt, final, graph);
        }
    }
}