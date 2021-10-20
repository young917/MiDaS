#include "setgenerator.hpp"

SetGenerator::SetGenerator(HyperGraph *graph) {
    this->node_core.resize(graph->number_of_nodes, -1);
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
        bin[vdeg]++;
    }
    int start = 0;
    for (int d = 0 ; d <= md ; d++){
        int num = bin[d];
        bin[d] += start;
        start += num;
    }
    pos.resize(graph->number_of_nodes, -1);
    vert.resize(graph->number_of_nodes, -1);
    for (int v = 0 ; v < graph->number_of_nodes ; v++){
        pos[v] = bin[node_core[v]];
        vert[pos[v]] = v;
        bin[node_core[v]]++;
    }
    for (int d = 1 ; d <= md ; d++){
        bin[d] = bin[d-1];
    }
    bin[0] = 0;
    
    check_hyperedge.resize(graph->number_of_hedges, false);
    cout << "Number of Nodes = " << graph->number_of_nodes << endl;
    for (int i = 0 ; i < graph->number_of_nodes ; i++){
        int v = vert[i];
        cout << "Set Node Core " << v;
        cout << " = " << node_core[v] << endl;
        int vdeg = (int)graph->node2hyperedge[v].size();
        for (int hidx = 0; hidx < vdeg ; hidx++){
            int h = graph->node2hyperedge[v][hidx];
            if (!check_hyperedge[h]){
                int hsize = (int)graph->hyperedge2node[h].size();
                for (int nvidx = 0 ; nvidx < hsize ; nvidx++){
                    int nv = graph->hyperedge2node[h][nvidx];
                    if (node_core[nv] > node_core[v]){
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
}

void SetGenerator::write_file(string outputdir, string dataname){
    string writeFile = outputdir + dataname + "_nodecoreness.txt";
    ofstream resultFile(writeFile.c_str());
    for (int v = 0 ; v < (int)node_core.size() ; v++){
        resultFile << v << "," << node_core[v] << endl;
    }
    resultFile.close();
}