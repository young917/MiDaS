from collections import defaultdict
from numpy.linalg import svd
from itertools import chain
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.linalg import svds, eigs, norm
import argparse
from tqdm import tqdm
import numpy as np
from numpy import linalg as LA
import math
import os
import matplotlib.pyplot as plt
import copy
plt.rcParams.update({'font.size': 12})

exceptdatas = ["threads-ask-ubuntu", "coauth-MAG-Geology-full", "coauth-MAG-History-full", "threads-math-sx", "coauth-DBLP-full"]
num = 300

def find_svs(dataname, inputpath, outputdir, temp_outputpath, recalculate_flag):
    outputpath = outputdir + "singular_values_full.txt"
    node2edges = defaultdict(list)
    nodename2index = {}
    node_index = 0
    number_of_nodes = 0
    number_of_edges = 0

    with open(inputpath, "r") as f:
        for idx, line in enumerate(f.readlines()):
            line = line[:-1] # strip enter
            nodes = line.split(",")
            for i in range(len(nodes)):
                node = nodes[i]
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
        return -1, -1, -1

    dim = min(number_of_edges, number_of_nodes)
    s = -1
    flag = True
    if (recalculate_flag is False) and (os.path.isfile(outputpath)):
        s = []
        with open(outputpath, "r") as f:
            for line in f.readlines():
                line = line[:-1]
                s.append(float(line))
        if (dataname in exceptdatas) is False and len(s) < dim:
            flag = True
        else:
            flag = False
    elif (recalculate_flag is False) and (dataname in exceptdatas):
        print("no singular value file", outputpath)
        flag = False
    elif len(temp_outputpath) > 0 and os.path.isfile(temp_outputpath):
        s = []
        with open(temp_outputpath, "r") as f:
            for line in f.readlines():
                line = line[:-1]
                order, singular_val = line.split(" ")
                s.append(float(singular_val))
        flag = False
    
    if (flag):
        try:
            incident_matrix = np.zeros(shape=(number_of_nodes, number_of_edges), dtype=np.byte)
            for v in node2edges.keys():
                for edge_idx in node2edges[v]:
                    incident_matrix[v,edge_idx] = 1
            # sum_of_squares = LA.norm(incident_matrix, 'fro') ** 2            
            s = svd(incident_matrix, compute_uv=False)
            s = sorted(list(s), reverse=True)
            assert len(s) == dim
            for i in range(1,dim):
                assert s[i-1] >= s[i]
            # assert math.fabs(sum_of_squares - np.sum(np.array(s) ** 2)) < 0.000001
            with open(outputpath, "w") as f:
                for _sv in s:
                    f.write(str(_sv) + "\n")
        except:
            # print("[" + inputpath + "] Error #V=" + str(number_of_nodes) + ", #E=" + str(number_of_edges))
            rows, cols = zip(*chain.from_iterable([[(v, edge_idx) for edge_idx in node2edges[v]] for v in node2edges.keys()]))
            nnz = len(rows)
            incident_matrix = coo_matrix((np.ones(nnz), (rows, cols)), shape=(number_of_nodes, number_of_edges))
            sum_of_squares = norm(incident_matrix, 'fro') ** 2
            rank = min(num , dim - 1)
            _, s, _ = svds(incident_matrix.tocsc(), k=rank)
            last_sv_square = sum_of_squares - sum([_s * _s for _s in s])
            last_sv = math.sqrt(last_sv_square)
            s = list(s)
            s.append(last_sv)
            assert len(s) == dim
            s = sorted(s, reverse=True)
            for i in range(1,dim):
                assert s[i-1] >= s[i]
            with open(outputpath, "w") as f:
                for _sv in s:
                    f.write(str(_sv) + "\n")

    with open(outputdir + "singular_values_info.txt", "w") as f:
        f.write("Dimension : " + str(dim) + "\n")
        # f.write("Sum of Squares : " + str(sum_of_squares) + '\n')

    return s, dim

def calculate_dstat(sampled_dict, answer_dict):
    # D-Statistics
    # Since svdist2 is not complete ... need min_max_key
    min_max_key = min( max(answer_dict.keys()) , max(sampled_dict.keys()) )
    keyset = set()
    keyset.update(list(answer_dict.keys()))
    keyset.update(list(sampled_dict.keys()))
    keys = list(keyset)
    keys = sorted(keys)

    answer_cumulsum = 0.0
    sampled_cumulsum = 0.0
    stat = 0.0
    for k in keys:
        if min_max_key < k:
            break
        if k in answer_dict:
            assert answer_cumulsum <= answer_dict[k]
            answer_cumulsum = answer_dict[k]
#             assert answer_cumulsum <= 1.0
        if k in sampled_dict:
            assert sampled_cumulsum <= sampled_dict[k]
            sampled_cumulsum = sampled_dict[k]
#             assert sampled_cumulsum <= 1.0, str(sampled_cumulsum)
        if stat < math.fabs(answer_cumulsum - sampled_cumulsum):
            stat = math.fabs(answer_cumulsum - sampled_cumulsum)
    return stat

def sv_dist2(_list_sv, dim, answer_max_portion):
    list_sv = copy.deepcopy(_list_sv)
    singular_values = {}
    
    if answer_max_portion != -1:
        number_of_required_svs = math.ceil(dim * answer_max_portion)
        list_sv = list_sv[:number_of_required_svs]
    denom = sum([_s * _s for _s in list_sv])
    until = 0.0
    total = dim
    for idx in range(len(list_sv)):
        proportion = (idx + 1) / total
        until += list_sv[idx] ** 2
        singular_values[proportion] = until / denom
#         if idx == (len(list_sv) - 1):
#             assert proportion == 1, proportion
#             assert singular_values[proportion] == 1, singular_values[proportion] 

    return singular_values

def show(answer_dict, sampled_dict, sampled_dir, stat, opt):
    plt.figure(dpi=120)
    plt.plot(answer_dict.keys(), answer_dict.values(), label="answer", linewidth=10, alpha=0.7)
    plt.plot(sampled_dict.keys(), sampled_dict.values(), label="sampled", linewidth=10, alpha=0.7)
    plt.legend()
    plt.title(str(stat))
    plt.tight_layout()
    plt.savefig(sampled_dir + "SVdist_" + str(opt) + ".jpg")
    plt.close()

if __name__ == "__main__":
    dataset = ["email-Enron-full", "email-Eu-full",
                "contact-primary-school", "contact-high-school",
                "NDC-substances-full", "NDC-classes-full", 
                "tags-math-sx", "tags-ask-ubuntu", "threads-ask-ubuntu", "coauth-MAG-Geology-full", "coauth-MAG-History-full"]
                #, "threads-math-sx", "coauth-DBLP-full"]
                
    sampled_algo_list = ["es_grid_search", "es_skewness_search", 
                    "es_grid_search_ext", "es_skewness_search_ext",
                    "es_global_deg_avg_grid_search", "es_global_deg_max_grid_search",
                    "ns_grid_search", "ns/global_deg_1.0000", "es/global_deg_min_0.0000", "ns/global_deg_0.0000",
                    "tihs", "rw/rw_c_1", "ff/ff_c_0.51_0.20",
                    "mh/add_degree", "mh/exchange_degree", "mh/remove_degree",
                    "mh/add_avg", "mh/exchange_avg", "mh/remove_avg"]
#                        "greedy/exchange/degree", "greedy/exchange/avg", "adjust"]

    parser = argparse.ArgumentParser()
    parser.add_argument('--repeat_str', required=False, default="1,2,3")
    parser.add_argument('--portion_str', required=False, default="0.30,0.20,0.10,0.40,0.50")
    parser.add_argument('--range_str', required=False, default="-1,-1")
    parser.add_argument('--dataname', required=False, default="")
    parser.add_argument('--algoname', required=False, default="")
    parser.add_argument('--recalculate', required=False, action='store_true')
    parser.add_argument('--use_temp', action='store_true')
    
    args = parser.parse_args()
    repeat_list = args.repeat_str.split(",")
    portion_list = args.portion_str.split(",")
    range_start, range_end = args.range_str.split(",")
    range_start, range_end = int(range_start), int(range_end)
    if len(args.dataname) > 0:
        dataset = args.dataname.split(",")
    if len(args.algoname) > 0:
        sampled_algo_list = args.algoname.split(",")
    if range_start != -1:
        dataset = dataset[range_start:range_end]

    for dataname in dataset:
        print("[[ " + dataname + " ]]")
        # for dataname in dataset:
        inputpath = "../dataset/" + dataname + ".txt"
        outputdir = "./results/answer_dist/" + dataname + "/"
        temp_outputpath = ""
        if args.use_temp:
            temp_outputpath = "./results/answer_dist/" + dataname + "/singular_values.txt"

        answer_s, answer_dim = find_svs(dataname, inputpath, outputdir, temp_outputpath, args.recalculate)
        if answer_s == -1:
            continue
        # answer_dict1 = sv_dist1(answer_s)
        answer_dict2 = sv_dist2(answer_s, answer_dim, -1)
        answer_max_portion = max(answer_dict2.keys())
        # answer_rep = representative_rate(answer_s, answer_sum_of_squares, answer_dim)
        print(answer_dim)
        with open("./results/answer_dist/" + dataname + "/singular_value_dist2.txt", "w") as f:
            for k,v in answer_dict2.items():
                f.write(str(k) + " " + str(v) + "\n")
        # with open("./results/answer_dist/" + dataname + "/singular_value_dist3.txt", "w") as f:
        #     f.write(str(answer_rep) + "\n")

        for portion_str in portion_list:
            for sampled_algo in tqdm(sampled_algo_list, desc=dataname + " " + portion_str):
                 for repeat_str in repeat_list:
                    # for dataname in dataset:
                    if "greedy/exchange" in sampled_algo:
                        sampled_dir = "./results/" + sampled_algo + "/" + dataname + "_" + portion_str + "/"
                    else:
                        sampled_dir = "./results/" + sampled_algo + "/" + dataname + "_" + portion_str + "/" + repeat_str + "/"
                    inputpath = sampled_dir + "sampled_graph.txt"
                    outputdir = sampled_dir
                    temp_outputpath = ""
                    if args.use_temp:
                        temp_outputpath = sampled_dir + "singular_values.txt"

                    if os.path.isfile(inputpath) is False:
                        with open("not_exist.txt", "a+") as f:
                            f.write(inputpath + "\n")
                        continue
                    elif dataname in exceptdatas and os.path.isfile(outputdir + "singular_values_full.txt") is False:
                        continue
                    else:
                        sampled_s, sampled_dim = find_svs(dataname, inputpath, outputdir, temp_outputpath, args.recalculate)
                    
                    sampled_dict2 = sv_dist2(sampled_s, sampled_dim, answer_max_portion)
                    with open(sampled_dir + "singular_value_dist2.txt", "w") as f:
                        for k,v in sampled_dict2.items():
                            f.write(str(k) + " " + str(v) + "\n")
                    
                    stat2 = calculate_dstat(sampled_dict2, answer_dict2)
                    with open(sampled_dir + "singular_value2_eval.txt", "w") as f:
                        f.write(str(stat2) + "\n")
                    
                    if "greedy/exchange" in sampled_algo:
                        break
