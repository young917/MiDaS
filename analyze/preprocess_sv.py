from collections import defaultdict
from tqdm import tqdm
import math
import os
import shutil
import argparse

def find_incidence_mat(inputpath, outputname):
    node2edges = defaultdict(list)
    nodename2index = {}
    node_index = 0
    number_of_nodes = 0
    number_of_edges = 0
    if os.path.isfile(inputpath) is False:
        with open("not_exist.txt", "+a") as f:
            f.write(inputpath + "\n")
        return -1
         
    with open(inputpath, "r") as f:
        for idx, line in enumerate(f.readlines()):
            line = line[:-1] # strip enter
            nodes = line.split(",")
            for i, node in enumerate(nodes):
                if node not in nodename2index:
                    nodename2index[node] = node_index
                    node_index += 1
                nodes[i] = nodename2index[node]
            for v in nodes:
                node2edges[int(v)].append(idx)
            number_of_edges += 1
    number_of_nodes = len(node2edges.keys())

    if number_of_edges == 0:
        print("Empty Input")
        return -1

    print(number_of_nodes, number_of_edges)
    total = 0
    with open(outputname + "dim.txt", "w") as f:
        f.write(str(number_of_nodes) + "\n")
        f.write(str(number_of_edges) + "\n")
    with open(outputname + "icd_row.txt", "w") as fr, open(outputname + "icd_col.txt", "w") as fc:
        for v in tqdm(range(number_of_nodes)):
            for h in node2edges[v]:
                fr.write(str(v + 1) + "\n")
                fc.write(str(h + 1) + "\n")
                total += 1
    with open(outputname + "icd.txt", "w") as f:
        for v in tqdm(range(number_of_nodes)):
            for h in node2edges[v]:
                f.write(str(v+1) + ":" + str(h+1) + "\n")
                total += 1
    print(total)
            
if __name__ == "__main__":
    dataset = ["threads-ask-ubuntu", "coauth-MAG-Geology-full", "coauth-MAG-History-full"]
    sampled_algo_list = ["midas", "ns/global_deg_1.0000", 
                    "es/global_deg_min_0.0000", "ns/global_deg_0.0000",
                    "tihs", "rw/rw_c_1", "ff/ff_c_0.51_0.20",
                    "hrw/noback", "hrw/skip",
                    "mgs/add_degree", "mgs/exchange_degree", "mgs/remove_degree",
                    "mgs/add_avg", "mgs/exchange_avg", "mgs/remove_avg",
                    "midas_grid_ablation", "maxdegree", "avgdegree", "midas_ns"]


    parser = argparse.ArgumentParser()
    parser.add_argument('--repeat_str', required=False, default="1,2,3")
    parser.add_argument('--portion_str', required=False, default="0.30,0.20,0.10,0.40,0.50")
    parser.add_argument('--dataname', required=False, default="")
    parser.add_argument('--algoname', required=False, default="")
    parser.add_argument('--recalculate', action='store_true')

    
    args = parser.parse_args()
    repeat_list = args.repeat_str.split(",")
    portion_list = args.portion_str.split(",")
    if len(args.dataname) > 0:
        dataset = args.dataname.split(",")
    if len(args.algoname) > 0:
        sampled_algo_list = args.algoname.split(",")
    
    for dataname in dataset:
        inputpath = "../dataset/" + dataname + ".txt"
        outputdir = "input/answer_dist/" + dataname + "/"
        if os.path.isdir(outputdir) is False:
            os.makedirs(outputdir)
        
        if (args.recalculate is False) and (os.path.isfile("output/answer_dist/" + dataname + "/singular_values_full.txt")):
            #(os.path.isfile("../Hypergraph_Sampling_cpp/results/answer_dist/" + dataname + "/singular_values_full.txt")):
            shutil.rmtree(outputdir)
        else:    
            later_outputdir = "output/answer_dist/" + dataname + "/"
            if os.path.isdir(later_outputdir) is False:
                os.makedirs(later_outputdir)
            answer_s = find_incidence_mat(inputpath, outputdir)

        for portion_str in portion_list:
            for algoname in tqdm(sampled_algo_list, desc=dataname + " " + portion_str):        
                for repeat_str in repeat_list:
                    inputpath = "../Hypergraph_Sampling_cpp/results/" + algoname + "/" + dataname + "_" + portion_str + "/" + repeat_str + "/sampled_graph.txt"
                    outputdir = "input/" + algoname + "/" + dataname + "_" + portion_str + "/" + repeat_str + "/"
                    if os.path.isdir(outputdir) is False:
                            os.makedirs(outputdir)
                    if (args.recalculate is False) and (os.path.isfile("output/" + algoname + "/" + dataname + "_" + portion_str + "/" + repeat_str + "/singular_values_full.txt")):
                        #(os.path.isfile("../Hypergraph_Sampling_cpp/results/" + algoname + "/" + dataname + "_" + portion_str  + "/" + repeat_str + "/singular_values_full.txt")):
                        shutil.rmtree(outputdir)
                    else:    
                        later_outputdir = "output/" + algoname + "/" + dataname + "_" + portion_str + "/" + repeat_str + "/"
                        if os.path.isdir(later_outputdir) is False:
                            os.makedirs(later_outputdir)
                        answer_s = find_incidence_mat(inputpath, outputdir)
    
    # for dataname in dataset:
    #     # for dataname in dataset:
    #     inputpath = "../dataset/" + dataname + ".txt"
    #     answer_s = find_incidence_mat(inputpath)
