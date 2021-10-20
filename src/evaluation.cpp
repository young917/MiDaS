#include "evaluation.hpp"

Evaluation::Evaluation(vector<int> &hypergraph_masking, HyperGraph *graph, string outputdir){
    this->outputdir = outputdir;
    int vindex = 0;
    int hindex = 0;
    map<int, int> node2index;
    
    int count = 0;
    for(int h = 0 ; h < graph->number_of_hedges ; h++){
        if(hypergraph_masking[h] == 1){
            count++;
        }
    }
    if (count == graph->number_of_hedges){
        hindex = graph->number_of_hedges;
        vindex = graph->number_of_nodes;
        node2hyperedge = graph->node2hyperedge;
        hyperedge2node = graph->hyperedge2node;
    }
    else{
        for(int h = 0 ; h < graph->number_of_hedges ; h++){
            if(hypergraph_masking[h] == 1){
                //cout << "[ " << h << " ]" << endl;
                int hsize = (int)graph->hyperedge2node[h].size();
                vector<int> tmp;
                for(int i = 0 ; i < hsize ; i++){
                    int v = graph->hyperedge2node[h][i];
                    if(node2index.find(v) == node2index.end()){
                        //cout << "add node2index" << endl;
                        node2index[v] = vindex++;
                        this->node2hyperedge.resize(vindex);
                    }
                    //cout << "add node2hyperedge" << endl;
                    tmp.push_back(node2index[v]);
                    this->node2hyperedge[node2index[v]].push_back(hindex);
                }
                //cout << "add hyperedge2node" << endl;
                this->hyperedge2node.push_back(tmp);
                hindex++;
            }    
        }
    }
    this->number_of_hedges = hindex;
    this->number_of_nodes= vindex;

    cout << "Number of Hyperedges : " << hindex << endl;
    cout << "Number of Nodes : " << vindex << endl;

    // clique-expansion
    if (this->number_of_nodes < 6000){
        this->node2node.resize(vindex);
        for(int va = 0 ; va < vindex ; va++){
            for(auto h : this->node2hyperedge[va]){
                int hsize = (int)this->hyperedge2node[h].size();
                for(int i = 0 ;  i < hsize; i++){
                    int vb = this->hyperedge2node[h][i];
                    if(va != vb){
                        this->node2node[va][vb] = 1;
                        this->node2node[vb][va] = 1;
                    }
                }
            }
        }
    }
}

void Evaluation::get_clustering_coef(){
    string writeFile = outputdir + "clusteringcoef.txt";
    cout << "Output Clustering Coefficient " << outputdir << endl;
    ofstream resultFile(writeFile.c_str());
    vector<int> neighbors;
    double cc = 0.0;
    double denominator = 0.0;
    for(int v = 0 ;  v < number_of_nodes; v++){
        cc = 0.0;
        denominator = 0.0;
        neighbors.clear();
        for(int nv = 0 ;  nv < number_of_nodes ; nv++){
            if(node2node[v][nv] == 1){
                neighbors.push_back(nv);
            }
        }
        int neighbor_size = (int)neighbors.size();
        for(int va_idx = 0 ; va_idx < neighbor_size ; va_idx++){
            for(int vb_idx = va_idx + 1 ; vb_idx < neighbor_size ; vb_idx++){
                int va = neighbors[va_idx];
                int vb = neighbors[vb_idx];
                cc += node2node[va][vb];
                denominator++;
            }
        }
        if (denominator != 0){
            cc /= denominator;
            resultFile << to_string(cc) << endl; // output for each node
        }
    }
    resultFile.close();
}

void Evaluation::get_global_clustering_coef(){
    random_device rd;
    mt19937 g(rd());
    
    int closed_wedges = 0;
    int wedges = 0;
    for (int t = 0; t < 1000000 ; t++){
        int rand_v = rand() % number_of_nodes;
        set<int> neighbors;
        int vdeg = (int)node2hyperedge[rand_v].size();
        for (int hidx = 0; hidx < vdeg ; hidx++){
            int h = node2hyperedge[rand_v][hidx];
            int hsize = (int)hyperedge2node[h].size();
            for (int nvidx = 0; nvidx < hsize ; nvidx++){
                int nv = hyperedge2node[h][nvidx];
                if (rand_v != nv){
                    neighbors.insert(nv);
                }
            }
        }
        if ((int)neighbors.size() < 2){
            continue;
        }
        vector<int> neighbors_vec;
        for (auto nv : neighbors){
            neighbors_vec.push_back(nv);
        }
        shuffle(neighbors_vec.begin(), neighbors_vec.end(), default_random_engine(rd()));
        
        int rand_nv1 = neighbors_vec[0];
        int rand_nv2 = neighbors_vec[1];
        bool find = false;
        vdeg = (int)node2hyperedge[rand_nv1].size();
        for (int hidx = 0; hidx < vdeg ; hidx++){
            int h = node2hyperedge[rand_nv1][hidx];
            int hsize = (int)hyperedge2node[h].size();
            for (int nvidx = 0; nvidx < hsize ; nvidx++){
                int nv = hyperedge2node[h][nvidx];
                if (nv == rand_nv2){
                    find = true;
                    closed_wedges++;
                    break;
                }
            }
            if (find){
                break;
            }
        }
        wedges++;
    }
    double global_cc = (double)closed_wedges / wedges;
    string writeFile = outputdir + "global_cc.txt";
    cout << "Output Global Clustering Coefficient " << outputdir << endl;

    ofstream resultFile(writeFile.c_str());
    resultFile << to_string(global_cc) << endl;
    resultFile.close();

    string compareFileName = outputdir + "clusteringcoef.txt";
    if (file_exist(compareFileName)){
        ifstream compareFile(compareFileName.c_str());
        double sum_cc = 0.0;
        int num = 0;
        string line;
        while (getline(compareFile, line)){
            double cc = stof(line); 
            sum_cc += cc;
            num++;
        }
        double avg_cc = sum_cc / num; 
        compareFile.close();
        
        cout << "Global CC ; Avg CC = " << to_string(global_cc) << " ; " << to_string(avg_cc) << endl;
    }
}

void Evaluation::get_effective_diameter(string path){
    // prepare "tmp.txt"'
    string writeFile = "./temporary/" + path + ".txt";
    ofstream resultFile(writeFile.c_str(), ofstream::out);
    resultFile << outputdir << endl;
    cout << "Make temporary file to get effective diameter" << endl;
    for (int h = 0 ; h < number_of_hedges ; h++){
        int hsize= (int)hyperedge2node[h].size();
        for (int i = 0 ;  i < hsize ; i++){
            int v = hyperedge2node[h][i];
            resultFile << v;
            if (i != (hsize - 1)){
                resultFile << ",";
            }
        }
        resultFile << endl;
    }
    resultFile.close();

    /*
    # include "Python.h"
    PyObject* get_diam = PyImport_ImportModule("helper.get_diameter");
    if(get_diam) {
        PyObject* diam = PyObject_GetAttrString(get_diam, "get_diameter");
        if(diam) {
            cout << "Call Python code" << endl;
            PyObject* r = PyObject_CallFunction(diam, NULL);
            if (r){
                int result = (int)PyInt_AS_LONG(r);
                Py_XDECREF(r);
            }
            Py_XDECREF(diam);
        }
        Py_XDECREF(get_diam);
    }

    string writeFile2 = outputdir + "eff_diameter.txt";
    ofstream resultFile2(writeFile2.c_str(), "w");
    resultFile2 << result << endl;
    resultFile2.close();
    */
}

void Evaluation::get_path_length(){
    vector<vector<int>> dist(number_of_nodes, vector<int>(number_of_nodes, number_of_nodes + 1));
    for(int v = 0 ; v < number_of_nodes ; v++){
        for(int nv = v + 1 ; nv < number_of_nodes ; nv++){
            if(node2node[v][nv] == 1){
                dist[v][nv] = node2node[v][nv];
                dist[nv][v] = node2node[v][nv];
            }
        }
        dist[v][v] = 0;
    }
    for(int k = 0 ; k < number_of_nodes ; k++){
        for(int i = 0 ; i < number_of_nodes ; i++){
            for(int j = 0 ; j < number_of_nodes ; j++){
                if(dist[i][j] > (dist[i][k] + dist[k][j])){
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }
    map<int, long long> pathdist;
    for(int v = 0 ; v < number_of_nodes ; v++){
        for(int nv = v + 1 ; nv < number_of_nodes ; nv++){
            int pathlength = dist[v][nv];
            pathdist[pathlength]++;
        }
    }
    string writeFile = outputdir + "pathlength.txt";
    cout << "Output Path Length " << outputdir << endl;
    ofstream resultFile(writeFile.c_str());
    resultFile << "pathlength,frequency" << endl;
    for(auto p : pathdist){
        if(p.first == (number_of_nodes + 1)){
            continue;    
        }
        resultFile << p.first << "," << p.second << endl;
    }
    resultFile.close();
}

int Evaluation::bfs(vector<bool> &check, int start_node){
    int number_of_visited = 0;
    queue<int> visited;
    visited.push(start_node);
    check[start_node] = true;
    number_of_visited++;
    while((int)visited.size() > 0){
        int v = visited.front();
        visited.pop();
        int vdeg = (int)node2hyperedge[v].size();
        for (int hidx = 0; hidx < vdeg ; hidx++){
            int h = node2hyperedge[v][hidx];
            int hsize = (int)hyperedge2node[h].size();
            for (int nvidx = 0; nvidx < hsize ; nvidx++){
                int nv = hyperedge2node[h][nvidx];
                if (check[nv] == false){
                    visited.push(nv);
                    check[nv] = true;
                    number_of_visited++;
                }
            }
        }
    }
    return number_of_visited;
}

void Evaluation::count_wcc(){
    vector<int> wcc_sizes;
    vector<bool> check(number_of_nodes, false);
    for(int start = 0 ; start < number_of_nodes ; start++){
        if(check[start] == false){
            int size_wcc = bfs(check, start);
            wcc_sizes.push_back(size_wcc);
        }
    }
    sort(wcc_sizes.begin(), wcc_sizes.end(), greater<int>());

    string writeFile = outputdir + "sizewcc.txt";
    cout << "Output Size of Weakly Connected Component " << outputdir << endl;
    ofstream resultFile(writeFile.c_str());
    resultFile << "size_wcc,num_nodes" << endl;

    int cumul_nodes = 0;
    for(auto wcc_sz : wcc_sizes){
        cumul_nodes += wcc_sz;
        resultFile << to_string(cumul_nodes) << "," << to_string(number_of_nodes) << endl;
    }
    resultFile.close();
}

void Evaluation::get_edge_density(){
    string writeFile = outputdir + "density.txt";
    cout << "Output Density of graph " << outputdir << endl;
    ofstream resultFile(writeFile.c_str());
    resultFile << "num nodes,num edges"  << endl;
    resultFile << to_string(number_of_nodes) << "," << to_string(number_of_hedges) << endl;
    resultFile.close();    
}