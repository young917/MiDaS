# include <random>
# include <ctime>
# include "hset.hpp"

#ifndef TEST_HPP
#define TEST_HPP

using namespace std;

inline void test_dynamic_change(HyperGraph *graph, int num_tries){
    vector<int> hyperedges;
    for (int h = 0; h < graph->number_of_hedges ; h++){
        hyperedges.push_back(h);
    }
    random_device rd;
    mt19937 g(rd());
    set<int> tmp;
    srand(static_cast<unsigned int>(std::time(0)));
    
    printf("Add All\n"); 
    shuffle(hyperedges.begin(), hyperedges.end(), g);
    tmp.insert(hyperedges[0]);
    HSet *testing1 = new HSet(tmp, graph, "avg");
    for(int i = 1 ; i < graph->number_of_hedges ; i++){
        tmp.clear();
        tmp.insert(hyperedges[i]);
        testing1->dynamic_update_eval(tmp, graph, "+");
    }
    for (auto e : testing1->evaluation){
        if(e.second != 0.0){
            cout << e.first << " " << e.second << endl;
        }
        assert(e.second == 0.0);
    }

    printf("==================================================================\n");
    printf("Delete All\n");
    shuffle(hyperedges.begin(), hyperedges.end(), g);
    tmp.clear();
    for(int h = 0 ; h < graph->number_of_hedges ; h++) tmp.insert(h);
    HSet *testing2 = new HSet(tmp, graph, "avg");
    for(int i = 0 ; i < graph->number_of_hedges ; i++){
        tmp.clear();
        tmp.insert(hyperedges[i]);
        testing2->dynamic_update_eval(tmp, graph, "-");
    }
    printf("Result\n");
    for (auto e : testing2->evaluation){
        if ( (1.0 - e.second) >= 1e-6){
            printf("%s : %.5f    \n", e.first.c_str(), e.second);
        }
        assert( (1.0 - e.second) < 1e-6 );
    }
    printf("\n");

    printf("==================================================================\n");
    printf("Randomly Add or Delete\n");
    shuffle(hyperedges.begin(), hyperedges.end(), g);
    int tmp_size = rand() % graph->number_of_hedges;
    set<int> current_set;
    tmp.clear();
    for (int i = 0 ; i < tmp_size ; i++){
        tmp.insert(hyperedges[i]);
        current_set.insert(hyperedges[i]);
    }
    HSet *testing3 = new HSet(tmp, graph, "avg");
    
    for (int i = 0; i< num_tries; i++){
        int selection = rand() % 2;
        if (((int)current_set.size() == 0) || (selection == 0)){// Add
            vector<int> add_pool;
            for(int h = 0; h < graph->number_of_hedges ; h++){
                if(testing3->hypergraph_masking[h] == 0) add_pool.push_back(h);
            }
            int s = rand() % (int)add_pool.size();
            shuffle(add_pool.begin(), add_pool.end(), g);
            tmp.clear();
            for(int j = 0; j < s ; j++){
                tmp.insert(add_pool[j]);
                current_set.insert(add_pool[j]);
            }
            testing3->dynamic_update_eval(tmp, graph, "+");
        }
        else if (((int)current_set.size() == graph->number_of_hedges) || (selection == 1)){// Delete
            vector<int> delete_pool;
            for(int h = 0; h < graph->number_of_hedges ; h++){
                if(testing3->hypergraph_masking[h] == 1) delete_pool.push_back(h);
            }
            int s = rand() % (int)delete_pool.size();
            shuffle(delete_pool.begin(), delete_pool.end(), g);
            tmp.clear();
            for(int j = 0; j < s ; j++){
                tmp.insert(delete_pool[j]);
                current_set.erase(delete_pool[j]);
            }
            testing3->dynamic_update_eval(tmp, graph, "-");
        }
        HSet *ans = new HSet(current_set, graph, "avg");
        
        for (auto dist : testing3->attr->attr_list){
            string dist_name = dist.first;
            for (auto d: dist.second){
                int k = d.first;
                int v = d.second;
                if (v != ans->attr->attr_list[dist_name][k]) assert(false);
            }
            for (auto d: ans->attr->attr_list[dist_name]){
                int k = d.first;
                int v = d.second;
                if (dist.second[k] != v) assert(false);
            }
        }
    }
    printf("==================================================================\n");
    printf("\nSuccess\n");
}
#endif