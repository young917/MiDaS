# include <random>
# include <ctime>
# include <sys/types.h>
# include "hset.hpp"

#ifndef SEARCH_HPP
#define SEARCH_HPP

using namespace std;

inline void query(string dirname, HyperGraph *graph, double &Dstat, double &time){
    // Read result files in repetition, then check whether errors exist

    int count = 0;
    Dstat = 0.0;
    time = 0;

    for (int i = 1 ; i < 4 ; i++){
        if (!file_exist(dirname + "/" + to_string(i) + "/sampled_graph.txt")){
            cout << "No exist " + dirname + "/" + to_string(i) + "/sampled_graph.txt" << endl;
            continue;
        }
        // HSet* queryHSet = new HSet(dirname + "/" + to_string(i) + "/", graph, "degree");
        int queryms; 
        double querydstat;
        set<int> temp;
        HSet* helper = new HSet(temp, graph, "degree");
        
        cout << "Read Sampled Graph" << endl;
        vector< vector<int> > samplehyperedge2node;
        unordered_map<string, int> nodename2index;
        ifstream samplegraphFile(dirname + "/" + to_string(i) + "/sampled_graph.txt");
        string line;
        int num_nodes = 0;

        while (getline(samplegraphFile, line)){
            vector<string> nodes = split(line, ',');
            vector<int> tokens;
            for (int i = 0; i < (int)nodes.size(); i++){
                if (nodename2index.find(nodes[i]) == nodename2index.end()){
                    int index = num_nodes++;
                    nodename2index[nodes[i]] = index;
                }
                int node_index = nodename2index[nodes[i]]; // string -> integer
                tokens.push_back(node_index);
            }
            sort(tokens.begin(), tokens.end());
            samplehyperedge2node.push_back(tokens);
        }
        cout << "Calculate Dist" << endl;

        vector<int> node_degree(num_nodes, 0);
        unordered_map<int, long long> degree;
        
        auto start = std::chrono::steady_clock::now();
        // queryHSet->calculate_distribution(graph, "degree");
        // queryHSet->evaluate(graph->attr);
        for(int hidx = 0; hidx < (int)samplehyperedge2node.size(); hidx++){
            int hsize = (int)samplehyperedge2node[hidx].size();
            for(int i = 0 ; i < hsize ; i++){
                int v = samplehyperedge2node[hidx][i];
                node_degree[v]++;
            }
        }
        for (int v = 0; v < num_nodes ; v++){
            int d = node_degree[v];
            if (d != 0){
                degree[d]++;
            }
            assert (d <= (int)samplehyperedge2node.size());
        }
        querydstat = helper->get_Dstat(degree, graph->attr->attr_list["degree"], "degree");
        auto end = std::chrono::steady_clock::now();
        queryms = std::chrono::duration_cast<chrono::milliseconds>(end - start).count();
        cout << queryms << endl;

        cout << "Read Results" << endl;

        string resultpath = dirname +  "/" + to_string(i) + "/result.txt";
        ifstream resultFile(resultpath.c_str());
        double check_dstat = 0;
        while (getline(resultFile, line)){
            vector<string> tokens = split(line, ':');
            // cout << tokens[0] << "," << endl;
            if (tokens[0].compare("Time ") == 0){
                string t = tokens[1];
                t.erase(std::remove(t.begin(), t.end(), '\n'), t.end());
                t.erase(std::remove(t.begin(), t.end(), ' '), t.end());
                cout << tokens[1] << endl;
                queryms += stoi(t);
            }
            else if (tokens[0].compare("degree ") == 0){
                string d = tokens[1];
                d.erase(std::remove(d.begin(), d.end(), '\n'), d.end());
                d.erase(std::remove(d.begin(), d.end(), ' '), d.end());
                check_dstat = stod(d);
            }
        }
        cout << to_string(check_dstat) << "\t" << to_string(querydstat) << endl;
        assert ( fabs(check_dstat - querydstat) < 1e-6);
        time += queryms;
        
        cout << "Read Degree Dist" << endl;

        // check: read degree dist and compare with queryHSet
        string degree_path = dirname +  "/" + to_string(i) + "/degree_dist.txt";
        ifstream degFile(degree_path.c_str());
        getline(degFile, line);
        while(getline(degFile, line)){
            vector<string> tokens = split(line, ',');
            int deg = stoi(tokens[0]);
            double freq = stod(tokens[1]);
            // cout << to_string((double)queryHSet->attr->attr_list["degree"][deg] / num_nodes) << "\t" << to_string(freq) << endl;
            assert ( fabs(((double)degree[deg] / num_nodes) - freq) < 1e-6 );
        }
        Dstat += querydstat;
        count++;
    }
    Dstat /= count;
    time /= count;
}

inline void adjust_timespent(HyperGraph *graph, double target_portion){
    vector<string> search_dirnames;
    stringstream stream;
    string portion_str;
    stream << std::fixed << std::setprecision(2) << target_portion;
    portion_str = stream.str();

    // ES Sampling 
    for (double i = -3 ; i < 6.5 ; i += 0.5){
        double alpha = pow((double)2, i);
        string alpha_str;
        
        stream.str("");
        stream << std::fixed << std::setprecision(4) << alpha;
        alpha_str = stream.str();

        string dirname = "results/es/global_deg_" + alpha_str + "/" + graph->dataname + "_" + portion_str;
        search_dirnames.push_back(dirname);
    }
    search_dirnames.push_back("results/es/global_deg_min_0.0000/" + graph->dataname + "_" + portion_str);

    // NS Sampling
    for (double i = -1 ; i <= 6 ; i += 1){
        double alpha = pow((double)2, i);
        string alpha_str;
        
        stream.str("");
        stream << std::fixed << std::setprecision(4) << alpha;
        alpha_str = stream.str();

        string dirname = "results/ns/global_deg_" + alpha_str + "/" + graph->dataname + "_" + portion_str;
        search_dirnames.push_back(dirname);
    }
    search_dirnames.push_back("results/ns/global_deg_0.0000/" + graph->dataname + "_" + portion_str);
    
    for (auto dirname : search_dirnames){
        vector<string> tokens = split(dirname, '/');
        string timedir  = "results_time";
        for (int i = 1 ; i < (int)(tokens.size()) ; i++){
            timedir += "/" + tokens[i];
            mkdir(timedir.c_str(), 0776);
        }

        string timepath = timedir +  "/adjust_time.txt";

        // if (file_exist(timepath)){
        //     ifstream timeFile(timepath.c_str());
        //     string line;
        //     getline(timeFile, line);
        //     if (stod(line) > 0.0){
        //         continue;
        //     }
        // }
        cout << dirname << endl;
        double qdstat;
        double qtime;
        query(dirname, graph, qdstat, qtime);
        
        ofstream timeFile(timepath.c_str());
        timeFile << to_string(qtime) << endl;
        timeFile.close();
    }
}
#endif