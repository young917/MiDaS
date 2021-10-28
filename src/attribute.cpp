#include "attribute.hpp"

Attribute::Attribute(vector<vector<int>> node2hyperedge, vector<vector<int>> hyperedge2node, string dataname){
    // All attributes do not include zero
    vector<int> vec;
    vector<int>::iterator it;
    int number_of_nodes = (int)node2hyperedge.size();
    int number_of_hedges = (int)hyperedge2node.size();

    // degree
    cout << "Get Degree Distribution" << endl;
    unordered_map<int, long long> degree;
    string degree_name = "./results/answer/" + dataname + "_degree.txt";
    if( file_exist(degree_name) ){
        // read file
        ifstream degFile(degree_name.c_str());
        string line;
        while (getline(degFile, line)){
            vector<string> tokens = split(line, ',');
            assert( (int)tokens.size() == 2);
            int deg = stoi(tokens[0]);
            long long freq = stoll(tokens[1]);
            degree[deg] = freq;
        }
    }
    else{
        for (int v = 0; v < number_of_nodes ; v++){
            int deg = (int)node2hyperedge[v].size();
            if(deg == 0){
                continue;
            }
            degree[deg]++;
        }
        // output file
        ofstream resultFile(degree_name.c_str());
        for (auto item : degree){
            resultFile << to_string(item.first) << "," << to_string(item.second) << endl;
        }
        resultFile.close();
    }
    this->attr_list["degree"] = degree;
    degree.clear();
    
    // size
    cout << "Get Size Distribution" << endl;
    unordered_map<int, long long> size;
    string size_name = "./results/answer/" + dataname + "_size.txt";
    if( file_exist(size_name) ){
        // read file
        ifstream sizeFile(size_name.c_str());
        string line;
        while (getline(sizeFile, line)){
            vector<string> tokens = split(line, ',');
            assert( (int)tokens.size() == 2);
            int sz = stoi(tokens[0]);
            long long freq = stoll(tokens[1]);
            size[sz] = freq;
        }
    }
    else{
        for (int h = 0; h < number_of_hedges ; h++){
            int hsize = (int)hyperedge2node[h].size();
            size[hsize]++;
        }
        // output file
        ofstream resultFile(size_name.c_str());
        for (auto item : size){
            resultFile << to_string(item.first) << "," << to_string(item.second) << endl;
        }
        resultFile.close();
    }
    this->attr_list["size"] = size;
    size.clear();

    // For Pair
    unordered_map<string, long long> tmp;
    unordered_map<int, long long> intersect;
    // intersect
    cout << "Get Intersect Distribution" << endl;
    string its_name = "./results/answer/" + dataname + "_intersect.txt";
    if( file_exist(its_name) ){
        // read file
        ifstream itsFile(its_name.c_str());
        string line;
        while (getline(itsFile, line)){
            vector<string> tokens = split(line, ',');
            assert( (int)tokens.size() == 2);
            int its = stoi(tokens[0]);
            long long freq = stoll(tokens[1]);
            intersect[its] = freq;
        }
    }
    else{
        for (int ha = 0; ha < number_of_hedges ; ha++){
            int hsize = (int)hyperedge2node[ha].size();
            for (int i = 0 ; i < hsize ; i++){
                int v = hyperedge2node[ha][i];
                int deg = (int)node2hyperedge[v].size();
                for (int j = 0 ; j < deg; j++){
                    int hb = node2hyperedge[v][j];
                    if (ha < hb){
                        string key = to_string(ha) + "_" + to_string(hb);
                        tmp[key]++;
                    }
                }
            }
        }
        for (auto item : tmp){
            intersect[item.second]++;
        }
        tmp.clear();
        // output file
        ofstream resultFile(its_name.c_str());
        for (auto item : intersect){
            resultFile << to_string(item.first) << "," << to_string(item.second) << endl;
        }
        resultFile.close();
    }
    this->attr_list["intersect"] = intersect;
    intersect.clear();

    // pairdeg
    cout << "Get PairDegree Distribution" << endl;
    unordered_map<int, long long> pairdeg;
    string pd_name = "./results/answer/" + dataname + "_pairdeg.txt";
    if( file_exist(pd_name) ){
        // read file
        ifstream pdFile(pd_name.c_str());
        string line;
        while (getline(pdFile, line)){
            vector<string> tokens = split(line, ',');
            assert( (int)tokens.size() == 2);
            int pd = stoi(tokens[0]);
            long long freq = stoll(tokens[1]);
            pairdeg[pd] = freq;
        }
    }
    else{
        for(int h = 0 ; h < number_of_hedges ; h++){
            int hsize = (int)hyperedge2node[h].size();
            for(int i = 0 ; i < hsize ; i++){
                for(int j = (i + 1) ; j < hsize ; j++){
                    int va = hyperedge2node[h][i];
                    int vb = hyperedge2node[h][j];
                    string key;
                    if (va > vb){
                        key = to_string(vb) + "_" + to_string(va);
                    }
                    else{
                        key = to_string(va) + "_" + to_string(vb);
                    }
                    tmp[key]++;
                }
            }
        }
        for(auto elem : tmp){
            pairdeg[elem.second]++;
        }
        tmp.clear();
        // output file
        ofstream resultFile(pd_name.c_str());
        for (auto item : pairdeg){
            resultFile << to_string(item.first) << "," << to_string(item.second) << endl;
        }
        resultFile.close();
    }
    this->attr_list["pairdeg"] = pairdeg;
    pairdeg.clear();
}

Attribute::Attribute(vector<vector<int>> node2hyperedge, vector<vector<int>> hyperedge2node, set<int> hset, vector<int> &node_degree, string eval_opt){
    vector<int> _hset; // set -> vector
    vector<int> _nodeset;

    for (auto h: hset){
        _hset.push_back(h);
        for (auto v: hyperedge2node[h]){
            if (find(_nodeset.begin(), _nodeset.end(), v) == _nodeset.end()){
                _nodeset.push_back(v);
            }
        }
    }
    sort(_hset.begin(), _hset.end());

    // degree
    unordered_map<int, long long> degree;
    for(auto h : hset){
        int hsize = (int)hyperedge2node[h].size();
        for(int i = 0 ; i < hsize ; i++){
            int v = hyperedge2node[h][i];
            node_degree[v]++;
        }
    }
    for (auto v : _nodeset){
        int d = node_degree[v];
        if (d != 0){
            degree[d]++;
        }
        assert(d <= (int)hyperedge2node.size());
    }
    this->attr_list["degree"] = degree;
    degree.clear();
    if (eval_opt.compare("avg") == 0){
        // size
        unordered_map<int, long long> size;
        for (auto h : _hset){
            int hsize = (int)hyperedge2node[h].size();
            size[hsize]++;
        }
        this->attr_list["size"] = size;
        size.clear();

        // For Pair
        unordered_map<string, long long> tmp;
        
        // intersect
        unordered_map<int, long long> intersect;
        for (int ha_idx = 0; ha_idx < (int)_hset.size() ; ha_idx++){
            int ha = _hset[ha_idx];
            int hsize = (int)hyperedge2node[ha].size();
            for (int i = 0 ; i < hsize ; i++){
                int v = hyperedge2node[ha][i];
                int deg = (int)node2hyperedge[v].size();
                for (int j = 0 ; j < deg; j++){
                    int hb = node2hyperedge[v][j];
                    if ( (hset.find(hb) != hset.end()) && (ha < hb)){
                        string key = to_string(ha) + "_" + to_string(hb);
                        tmp[key]++;
                    }
                }
            }
        }
        for (auto item : tmp){
            intersect[item.second]++;
        }
        tmp.clear();
        this->attr_list["intersect"] = intersect;
        intersect.clear();

        // pairdeg
        unordered_map<int, long long> pairdeg;
        for(int h_idx = 0 ; h_idx < (int)_hset.size() ; h_idx++){
            int h = _hset[h_idx];
            int hsize = (int)hyperedge2node[h].size();
            for(int i = 0 ; i < hsize ; i++){
                for(int j = (i + 1) ; j < hsize ; j++){
                    int va = hyperedge2node[h][i];
                    int vb = hyperedge2node[h][j];
                    string key;
                    if (va > vb){
                        key = to_string(vb) + "_" + to_string(va);
                    }
                    else{
                        key = to_string(va) + "_" + to_string(vb);
                    }
                    tmp[key]++;
                }
            }
        }
        for(auto elem : tmp){
            pairdeg[elem.second]++;
        }
        tmp.clear();
        this->attr_list["pairdeg"] = pairdeg;
        pairdeg.clear();
    }
}

Attribute::Attribute(Attribute *attr){
    unordered_map<int, long long> degree;
    unordered_map<int, long long> size;
    unordered_map<int, long long> intersect;
    unordered_map<int, long long> pairdeg;

    for (auto dist : attr->attr_list){
        string dist_name = dist.first;
        if (dist_name.compare("degree") == 0){
            assert((int)attr->attr_list.size() > 0);
            for(auto d : attr->attr_list["degree"]){
                degree[d.first] = d.second;
            }
            this->attr_list["degree"] = degree;
            degree.clear();
        }
        else if (dist_name.compare("size") == 0){
            for(auto d : attr->attr_list["size"]){
                size[d.first] = d.second;
            }
            this->attr_list["size"] = size;
            size.clear();
        }
        else if (dist_name.compare("intersect") == 0){
            for(auto d : attr->attr_list["intersect"]){
                intersect[d.first] = d.second;
            }
            this->attr_list["intersect"] = intersect;
            intersect.clear();
        }
        else if (dist_name.compare("pairdeg") == 0){
            for(auto d : attr->attr_list["pairdeg"]){
                pairdeg[d.first] = d.second;
            }
            this->attr_list["pairdeg"] = pairdeg;
            pairdeg.clear();
        }
    }
}

void Attribute::set_attribute(vector<vector<int>> node2hyperedge, vector<vector<int>> hyperedge2node, vector<int> hset, vector<int> &node_degree, string eval_opt){
    vector<int> _nodeset;
    fill(node_degree.begin(), node_degree.end(), 0);

    for (auto h: hset){
        for (auto v: hyperedge2node[h]){
            if (find(_nodeset.begin(), _nodeset.end(), v) == _nodeset.end()){
                _nodeset.push_back(v);
            }
        }
    }

    for (auto dist : attr_list){
        auto &d = dist.second;
        d.clear();
    }
    
    unordered_map<int, long long> degree;
    for(auto h : hset){
        int hsize = (int)hyperedge2node[h].size();
        for(int i = 0 ; i < hsize ; i++){
            int v = hyperedge2node[h][i];
            node_degree[v]++;
        }
    }
    for (auto v : _nodeset){
        int d = node_degree[v];
        if (d != 0){
            degree[d]++;
        }
        assert(d <= (int)hyperedge2node.size());
    }
    this->attr_list["degree"] = degree;
    degree.clear();
    if (eval_opt.compare("avg") == 0){
        // size
        unordered_map<int, long long> size;
        for (auto h : hset){
            int hsize = (int)hyperedge2node[h].size();
            size[hsize]++;
        }
        this->attr_list["size"] = size;
        size.clear();

        // For Pair
        unordered_map<string, long long> tmp;

        // intersect
        unordered_map<int, long long> intersect;
        for (int ha_idx = 0; ha_idx < (int)hset.size() ; ha_idx++){
            int ha = hset[ha_idx];
            int hsize = (int)hyperedge2node[ha].size();
            for (int i = 0 ; i < hsize ; i++){
                int v = hyperedge2node[ha][i];
                int deg = (int)node2hyperedge[v].size();
                for (int j = 0 ; j < deg; j++){
                    int hb = node2hyperedge[v][j];
                    auto it = find(hset.begin(), hset.end(), hb);
                    if ( (it != hset.end()) && (ha < hb)){
                        string key = to_string(ha) + "_" + to_string(hb);
                        tmp[key]++;
                    }
                }
            }
        }
        for (auto item : tmp){
            intersect[item.second]++;
        }
        tmp.clear();
        this->attr_list["intersect"] = intersect;
        intersect.clear();

        // pairdeg
        unordered_map<int, long long> pairdeg;
        for(int h_idx = 0 ; h_idx < (int)hset.size() ; h_idx++){
            int h = hset[h_idx];
            int hsize = (int)hyperedge2node[h].size();
            for(int i = 0 ; i < hsize ; i++){
                for(int j = (i + 1) ; j < hsize ; j++){
                    int va = hyperedge2node[h][i];
                    int vb = hyperedge2node[h][j];
                    string key;
                    if (va > vb){
                        key = to_string(vb) + "_" + to_string(va);
                    }
                    else{
                        key = to_string(va) + "_" + to_string(vb);
                    }
                    tmp[key]++;
                }
            }
        }
        for(auto elem : tmp){
            pairdeg[elem.second]++;
        }
        tmp.clear();
        this->attr_list["pairdeg"] = pairdeg;
        pairdeg.clear();
    }
}

void Attribute::to_textfile(string output_dir){
    for (auto a: attr_list){
        string writeFile  = output_dir + a.first + "_dist.txt";
        ofstream resultFile(writeFile.c_str());
        resultFile << a.first << ",value" << endl;
        double sum_val = 0.0;
        for (auto d: a.second){
            sum_val += d.second;
        }
        for (auto d: a.second){
            resultFile << to_string(d.first) << "," << to_string(d.second / sum_val) << endl;
        }
        resultFile.close();
    }
}

bool Attribute::operator==(Attribute& attr_p) const{
    for (auto dist : this->attr_list){
        string dist_name = dist.first;
        for (auto d: dist.second){
            int k = d.first;
            int v = d.second;
            if (v != attr_p.attr_list[dist_name][k]) return false;
        }
        for (auto d: attr_p.attr_list[dist_name]){
            int k = d.first;
            int v = d.second;
            if (dist.second[k] != v) return false;
        }
    }
    return true;
}