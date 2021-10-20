#include "hset.hpp"

HSet::HSet(set<int> subhypergraph, HyperGraph *graph, string eval_opt) {
    this->hypergraph_masking.resize(graph->number_of_hedges, false);
    for (auto h : subhypergraph){
        this->hypergraph_masking[h] = true;
    }
    this->setsize = subhypergraph.size();
    this->node_degree.resize(graph->number_of_nodes);
    this->evaluation_value = -1.0; // initial value
    this->attr = new Attribute(graph->node2hyperedge, graph->hyperedge2node, subhypergraph, node_degree, eval_opt);
}

HSet::HSet(HSet *hset){
    // copy
    for(int i = 0 ; i < (int)hset->hypergraph_masking.size() ; i++){
        this->hypergraph_masking.push_back(hset->hypergraph_masking[i]);
    }
    this->setsize = hset->setsize;
    this->attr = new Attribute(hset->attr);
    this->evaluation_value = hset->evaluation_value;
    for(auto e: hset->evaluation){
        this->evaluation[e.first] = e.second;
    }
    int num_nodes = (int)hset->node_degree.size();
    this->node_degree.resize(num_nodes, 0);
    for (int v = 0 ; v < num_nodes ; v++){
        this->node_degree[v] = hset->node_degree[v];
    }
}
HSet::HSet(string dirpath, HyperGraph *graph, string eval_opt){
    set<int> subhypergraph;
    string hepath = dirpath + "sampled_hyperedges.txt";
    if (file_exist(hepath)){
        ifstream graphFile(hepath.c_str());
        string line;
        while (getline(graphFile, line)){
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            int hidx = stoi(line);
            subhypergraph.insert(hidx);
        }
        graphFile.close();
    }
    else{
        string path = dirpath + "sampled_graph.txt";
        ifstream graphFile(path.c_str());
        string line;
        int linenum = 0;
        cout << "Read Sampled Graph" << endl; 
        while (getline(graphFile, line)){
            vector<string> tokens = split(line, ',');
            vector<int> nodes;
            for (int i = 0 ; i < (int)tokens.size() ; i++){
                int node = graph->nodename2index[tokens[i]];
                nodes.push_back(node);
            }
            auto hidx_addr = find(graph->hyperedge2node.begin(), graph->hyperedge2node.end(), nodes);
            assert(hidx_addr != graph->hyperedge2node.end());
            int hidx = hidx_addr - graph->hyperedge2node.begin();
            subhypergraph.insert(hidx);
            linenum += 1;
        }
        bool flag = ( ((int)subhypergraph.size() > 0) && (linenum == (int)subhypergraph.size()) );
        if (!flag){
            cout << "Directory is " << dirpath << endl;
        } 
        assert ( ((int)subhypergraph.size() > 0) && (linenum == (int)subhypergraph.size()) );
        graphFile.close();

        string outpath = dirpath + "sampled_hyperedges.txt";
        ofstream outFile(outpath.c_str());
        for (auto hidx : subhypergraph){
            outFile << to_string(hidx) << endl;
        }
        outFile.close();
    }
    this->hypergraph_masking.resize(graph->number_of_hedges, false);
    for (auto h : subhypergraph){
        this->hypergraph_masking[h] = true;
    }
    this->setsize = (int)subhypergraph.size();
    this->node_degree.resize(graph->number_of_nodes);
    this->evaluation_value = -1.0; // initial value
    this->attr = new Attribute(graph->node2hyperedge, graph->hyperedge2node, subhypergraph, node_degree, eval_opt);
}

vector<int> HSet::get_hyperedgeset(void){
    vector<int> hyperedgeset;
    for (int h = 0 ; h < (int)hypergraph_masking.size() ; h++){
        if (hypergraph_masking[h]) hyperedgeset.push_back(h);
    }
    sort(hyperedgeset.begin(), hyperedgeset.end());
    return hyperedgeset;
}

double HSet::get_Dstat(unordered_map<int,long long> &dist, unordered_map<int,long long> &dist_target, string eval_name){
    vector<int> key_list;
    vector<int>::iterator it;
    double cumul_dist_v = 0.0, cumul_dist_t_v = 0.0;
    double sum_dist = 0.0, sum_dist_t = 0.0;
    
    for (auto d : dist){
        if (d.second == 0) continue;
        sum_dist += d.second;
        it = find(key_list.begin(), key_list.end(), d.first);
        if (it == key_list.end()){
            key_list.push_back(d.first);
        }
    }
    for (auto d : dist_target){
        if (d.second == 0) continue;
        sum_dist_t += d.second;
        it = find(key_list.begin(), key_list.end(), d.first);
        if (it == key_list.end()){
            key_list.push_back(d.first);
        }
    }
    sort(key_list.begin(), key_list.end());

    double d_stat = 0.0;
    for (auto key : key_list){
        if (sum_dist != 0){
            cumul_dist_v += ((double)dist[key] / sum_dist);
        }
        if (sum_dist_t != 0){
            cumul_dist_t_v += ((double)dist_target[key] / sum_dist_t);
        }

        if (d_stat < fabs(cumul_dist_v - cumul_dist_t_v)){
            d_stat = fabs(cumul_dist_v - cumul_dist_t_v);
        }
    }
    return d_stat;
}

void HSet::calculate_distribution(HyperGraph *graph, string eval_opt){
    vector<int> hyperedgeset = get_hyperedgeset();
    attr->set_attribute(graph->node2hyperedge, graph->hyperedge2node, hyperedgeset, node_degree, eval_opt);
}

void HSet::evaluate(Attribute *target_attr){
    for (auto dist : attr->attr_list){
        string dist_name = dist.first;
        auto dist1 = dist.second;
        auto dist2 =  target_attr->attr_list[dist_name];
        evaluation[dist_name] = this->get_Dstat(dist1, dist2, dist_name);
    }
    double tmp = 0.0;
    for (auto eval : evaluation){
        tmp += eval.second;
    }
    evaluation_value = tmp / (double)evaluation.size();
}

void HSet::dynamic_degree(set<int> deltah, HyperGraph *graph, string sign){
    map<int, int> delta_deg;
    for (auto h : deltah){
        int hsize = (int)graph->hyperedge2node[h].size();
        for (int i = 0 ; i < hsize ; i++){
            int v = graph->hyperedge2node[h][i];
            delta_deg[v] += 1;
        }
    }
    for (auto delta : delta_deg){
        int v = delta.first;
        int d = delta.second;
        int prev_deg = node_degree[v];
        int cur_deg = 0;
        if (sign.compare("-") == 0) cur_deg = prev_deg - d;
        else if (sign.compare("+") == 0) cur_deg = prev_deg + d;
        node_degree[v] = cur_deg;

        if (prev_deg != 0) attr->attr_list["degree"][prev_deg] -= 1;
        if (cur_deg != 0) attr->attr_list["degree"][cur_deg] += 1;

        if ((prev_deg != 0) && (attr->attr_list["degree"][prev_deg] < 0)) cout << "Error Prev Deg: " << prev_deg << endl;
        if ((cur_deg != 0) && (attr->attr_list["degree"][cur_deg] < 0)) cout << "Error Cur Deg: " << cur_deg << endl;

        assert((prev_deg == 0) || (attr->attr_list["degree"][prev_deg] >= 0));
        assert((cur_deg == 0) || (attr->attr_list["degree"][cur_deg] >= 0));
    }
}

void HSet::dynamic_size(set<int> deltah, HyperGraph *graph, string sign){
    for (auto h : deltah){
        int hsize = (int)graph->hyperedge2node[h].size();
        if (sign.compare("-") == 0) attr->attr_list["size"][hsize] -= 1;
        else if (sign.compare("+") == 0) attr->attr_list["size"][hsize] += 1;
        assert(attr->attr_list["size"][hsize] >= 0);
    }
}

void HSet::dynamic_intersect(set<int> deltah, HyperGraph *graph, string sign){
    vector<bool> check(graph->number_of_hedges, false);
    vector<int> vec;
    vector<int>::iterator it;

    for (auto h : deltah){
        fill(check.begin(), check.end(), false);
        int hsize = (int)graph->hyperedge2node[h].size();
        for (int i = 0 ; i < hsize ; i++){
            int v = graph->hyperedge2node[h][i];
            int vdeg = (int)graph->node2hyperedge[v].size();
            for (int j = 0 ; j < vdeg ; j++){
                int nh = graph->node2hyperedge[v][j];
                if (h != nh && check[nh] == false){
                    check[nh] = true;
                    bool flag = false;
                    if ( (find(deltah.begin(), deltah.end(), nh) == deltah.end()) && (hypergraph_masking[nh] == 1) ) flag= true;
                    else if ( (find(deltah.begin(), deltah.end(), nh) != deltah.end()) && (h > nh) ) flag=true;
                    if (flag){
                        vec.resize(graph->number_of_hedges);
                        it = set_intersection(graph->hyperedge2node[h].begin(), graph->hyperedge2node[h].end(), graph->hyperedge2node[nh].begin(), graph->hyperedge2node[nh].end(), vec.begin());
                        int its_size = (int)(it - vec.begin());

                        if (its_size == 0) continue;
                        if (sign.compare("-") == 0) attr->attr_list["intersect"][its_size] -= 1;
                        else if (sign.compare("+") == 0) attr->attr_list["intersect"][its_size] += 1;
                        assert((its_size == 0) || (attr->attr_list["intersect"][its_size] >= 0));
                        vec.clear();
                    }
                }
            }
        }
    }
}

void HSet::dynamic_pairdegree(set<int> deltah, HyperGraph *graph, string sign){
    vector<int> before_hypergraph;
    vector<int> current_hyperedges = get_hyperedgeset();
    for (auto h : current_hyperedges) before_hypergraph.push_back(h);
    for (auto h : deltah){
        if (sign.compare("-") == 0) before_hypergraph.push_back(h);
        else if (sign.compare("+") == 0) before_hypergraph.erase(remove(before_hypergraph.begin(), before_hypergraph.end(), h), before_hypergraph.end());
    }
    sort(before_hypergraph.begin(), before_hypergraph.end());

    map<pair<int,int>, int> delta_deg;
    for (auto h : deltah){
        int hsize = (int)graph->hyperedge2node[h].size();
        for (int v_idx = 0 ; v_idx < hsize ; v_idx++){
            for (int nv_idx = v_idx + 1 ; nv_idx < hsize ; nv_idx++){
                int v = graph->hyperedge2node[h][v_idx];
                int nv = graph->hyperedge2node[h][nv_idx];
                pair<int, int> key;
                if (nv < v){
                    key.first = nv;
                    key.second = v;
                }
                else{
                    key.first = v;
                    key.second = nv;
                }
                delta_deg[key] += 1;
            }
        }
    }
    vector<int> vec;
    vector<int> vec2;
    vector<int>::iterator it;

    for (auto d : delta_deg){
        int v = d.first.first;
        int nv = d.first.second;
        int delta = d.second;
        vec.resize(graph->number_of_hedges);
        vec2.resize(graph->number_of_hedges);
        it = set_intersection(graph->node2hyperedge[v].begin(), graph->node2hyperedge[v].end(), graph->node2hyperedge[nv].begin(), graph->node2hyperedge[nv].end(), vec.begin());
        vec.resize(it - vec.begin());
        it = set_intersection(before_hypergraph.begin(), before_hypergraph.end(), vec.begin(), vec.end(), vec2.begin());
        int prev_deg = (int)(it - vec2.begin());
        int cur_deg = 0;
        if (sign.compare("-") == 0) cur_deg = prev_deg - delta;
        else if (sign.compare("+") == 0) cur_deg = prev_deg + delta;
        if (prev_deg != 0) attr->attr_list["pairdeg"][prev_deg] -= 1;
        if (cur_deg != 0) attr->attr_list["pairdeg"][cur_deg] += 1;
        assert((prev_deg == 0) || (attr->attr_list["pairdeg"][prev_deg] >= 0));
        assert((cur_deg == 0) || (attr->attr_list["pairdeg"][cur_deg] >= 0));
        vec.clear();
        vec2.clear();
    }
}

void HSet::dynamic_update_dist(set<int> deltaset, HyperGraph *graph, string sign){
    for (auto dist : attr->attr_list){
        string dist_name = dist.first;
        if (dist_name.compare("degree") == 0){
            dynamic_degree(deltaset, graph, sign);
        }
        else if (dist_name.compare("size") == 0){
            dynamic_size(deltaset, graph, sign);
        }
        else if (dist_name.compare("intersect") == 0){
            dynamic_intersect(deltaset, graph, sign);
        }
        else if (dist_name.compare("pairdeg") == 0){
            dynamic_pairdegree(deltaset, graph, sign);
        }
    }
}

void HSet::change_state(set<int> deltaset, string sign){
    for (auto h : deltaset){
        if (sign.compare("-") == 0){
            assert(hypergraph_masking[h] == 1);
            hypergraph_masking[h] = 0;
            setsize--;
        }
        else if (sign.compare("+") == 0){
            assert(hypergraph_masking[h] == 0);
            hypergraph_masking[h] = 1;
            setsize++;
        }
    }
}

double HSet::dynamic_update_eval(set<int> deltaset, HyperGraph *graph, string sign){
    // add/remove hyperedges & Update distribution dynamically & Get evaluation
    change_state(deltaset, sign);
    dynamic_update_dist(deltaset, graph, sign);
    evaluate(graph->attr);

    return evaluation_value;
}

double HSet::try_degree_update(set<int> deltaset, HyperGraph *graph, string sign){
    // Keep degree distribution in attr->attr_list["degree"]
    double ret = 0;
    change_state(deltaset, sign);
    dynamic_degree(deltaset, graph, sign);
    ret = get_Dstat(attr->attr_list["degree"], graph->attr->attr_list["degree"], "degree");
    
    // revert
    if (sign.compare("-") == 0) sign = "+";
    else if (sign.compare("+") == 0) sign = "-";
    change_state(deltaset, sign);
    dynamic_degree(deltaset, graph, sign);
    
    return ret;
}

double HSet::try_size_update(set<int> deltaset, HyperGraph *graph, string sign){
    // Keep degree distribution in attr->attr_list["degree"]
    double ret = 0;
    double prev_ret = get_Dstat(attr->attr_list["size"], graph->attr->attr_list["size"], "size");
    
    change_state(deltaset, sign);
    dynamic_size(deltaset, graph, sign);
    ret = get_Dstat(attr->attr_list["size"], graph->attr->attr_list["size"], "size");
    
    // revert
    if (sign.compare("-") == 0) sign = "+";
    else if (sign.compare("+") == 0) sign = "-";
    change_state(deltaset, sign);
    dynamic_size(deltaset, graph, sign);
    
    return prev_ret - ret;
}

bool HSet::operator==(HSet &hset_p) const{
    if (attr != hset_p.attr) return false;
    for (int h = 0 ; h < (int)hypergraph_masking.size() ; h++){
        if (hypergraph_masking[h] != hset_p.hypergraph_masking[h]) return false;
    }
    if (setsize != hset_p.setsize) return false;
    return true;
}

void HSet::save_as_txt(HyperGraph *graph, string outputdir, string outputname){
    set<int> nodeindexset;

    if (graph->exist_edgename){
        string sampledgraphCite = outputdir + outputname + ".cites";
        ofstream sampledgraphFile(sampledgraphCite.c_str());
        for(int h = 0 ; h < graph->number_of_hedges ; h++){
            if (hypergraph_masking[h] == 1){
                int hsize = (int)graph->hyperedge2node[h].size();
                for( int i = 0 ; i < hsize ; i++){
                    int v = graph->hyperedge2node[h][i];
                    nodeindexset.insert(v);
                    sampledgraphFile << graph->index2nodename[v] << "\t" << graph->edgename[h] << endl; 
                }
            }
        }
        string sampledgraphContent = outputdir + outputname + ".content_wtest";
        ofstream sampledgraphFile2(sampledgraphContent.c_str());
        string input_content = graph->inputpath + graph->dataname + "/" + graph->dataname + ".content_wtest"; 
        ifstream contentFile(input_content.c_str());
        string line;
        while (getline(contentFile, line)){
            vector<string> tmp = split(line, '\t');
            int nodeindex = graph->nodename2index[tmp[0]];
            if (nodeindexset.find(nodeindex) != nodeindexset.end()){
                // when this node is sampled
                sampledgraphFile2 << line << endl;
            }
        }
        sampledgraphFile2.close();
        contentFile.close();
    }
    else{
        string sampledgraphFName = outputdir + outputname + ".txt";
        ofstream sampledgraphFile(sampledgraphFName.c_str());
        for(int h = 0 ; h < graph->number_of_hedges ; h++){
            if (hypergraph_masking[h] == 1){
                int hsize = (int)graph->hyperedge2node[h].size();
                for( int i = 0 ; i < hsize ; i++){
                    int v = graph->hyperedge2node[h][i];
                    sampledgraphFile << graph->index2nodename[v]; 
                    if(i != (hsize -1)){
                        sampledgraphFile << ",";
                    }
                }
                sampledgraphFile << endl;
            }
        }
    }
}