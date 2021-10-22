import snap
import os
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.linalg import svds, eigs
import numpy as np
from itertools import chain
import random
from collections import defaultdict
import argparse
from tqdm import tqdm
import networkx as nx
import math
import time
import pandas as pd

dataset = ["email-Enron-full", "email-Eu-full",
            "contact-high-school", "contact-primary-school",
            "NDC-classes-full", "NDC-substances-full",
            "tags-ask-ubuntu", "tags-math-sx",
            "threads-ask-ubuntu", "coauth-MAG-History-full", "coauth-MAG-Geology-full"]
            # "threads-math-sx","coauth-DBLP-full"]
repeat_time = 3

def get_subdirectories(ls, dataname=None):
    ret = []
    cur = ls
    while True:
        ret += [a + "/" for a in cur if os.path.isfile(a + "/sampled_graph.txt")]
        next = [os.path.join(a,b) for a in cur for b in os.listdir(a) if os.path.isdir(os.path.join(a,b))]
        cur = next
        if len(next) == 0:
            break
    if dataname is not None:
        ret = [r for r in ret if dataname in r]
    return ret

def get_diameter(inputpath, dataname, algorithmlist, portionlist, split_num, split_idx, recalculate_flag):
    print("Analyzing Effective Diameter...")
    # read data
    if dataname is not None:
        _dataset = [dataname]
    else:
        _dataset = dataset

    if len(inputpath) > 0:
        # get nested all directories
        dir_list_all = get_subdirectories([inputpath])
    else:
        dir_list_all = []
        # answer
        for dataname in _dataset:
            inputpath = "../../dataset/" + dataname + ".txt"
            outputdir = "../results/answer_dist/" + dataname + "/"
            dir_list_all.append((inputpath, outputdir))
        # algorithms
        for algoname in algorithmlist:
            for dataname in _dataset:
                for portion in portionlist:
                    # repeat
                    for i in range(1, repeat_time + 1):
                        dir_name = "../results/" + algoname + "/" + dataname + "_" + portion + "/" + str(i) + "/" 
                        if os.path.isdir(dir_name):
                            dir_list_all.append(dir_name)
                        # for mhname in mhlist:
                        #     mhdir_name = "../results/" + algoname + "_" + mhname + "/" + dataname + "_" + portion + "/" + str(i) + "/"
                        #     if os.path.isdir(mhdir_name):
                        #         dir_list_all.append(mhdir_name)
                    # no repeat
                    dir_name = "../results/" + algoname + "/" + dataname + "_" + portion + "/"
                    if os.path.isdir(dir_name):
                        dir_list_all.append(dir_name)
                    # for mhname in mhlist:
                    #     mhdir_name = "../results/" + algoname + "_" + mhname + "/" + dataname + "_" + portion + "/"
                    #     if os.path.isdir(mhdir_name):
                    #         dir_list_all.append(mhdir_name)

        # for algoname in algorithmlist_nmh:
        #     for dataname in _dataset:
        #         for portion in portionlist:
        #             for i in range(1, repeat_time + 1):
        #                 dir_name = "../results/" + algoname + "/" + dataname + "_" + portion + "/" + str(i) + "/" 
        #                 if os.path.isdir(dir_name):
        #                     dir_list_all.append(dir_name)
        #             dir_name = "../results/" + algoname + "/" + dataname + "_" + portion + "/"
        #             if os.path.isdir(dir_name):
        #                 dir_list_all.append(dir_name)
    split_gap = math.ceil(len(dir_list_all) / split_num)
    end_idx = min(len(dir_list_all), split_gap * (split_idx + 1))
    dir_list_all = dir_list_all[split_gap * (split_idx) : end_idx]
    for dir_name in tqdm(dir_list_all):
        #print(dir_name)
        if type(dir_name) is tuple:
            inputpath, outputdir = dir_name
        else:
            inputpath = dir_name + "sampled_graph.txt"
            outputdir = dir_name

        flag_sampled_graph = os.path.isfile(inputpath)
        flag_effdiameter = os.path.isfile(outputdir + "effdiameter.txt")
        if flag_sampled_graph is False:
            continue
        elif (recalculate_flag is False) and (flag_effdiameter is True):
            continue

        edges = []
        nodename2index = {}
        node_index = 0
        with open(inputpath, "r") as f:
            for line in f.readlines():
                line = line[:-1]
                nodes = line.split(",")
                for i, node in enumerate(nodes):
                    if node not in nodename2index:
                        nodename2index[node] = node_index
                        node_index += 1
                    # reindex
                    nodes[i] = nodename2index[node]
                edges.append(nodes)

        if len(edges) == 0:
            print("Empty Input")
            return
        if os.path.isdir(outputdir) is False:
            continue

        # build snap graph
        pg = snap.TUNGraph.New()
        nodeset = set([])
        for hyperedge in edges:
            for n in hyperedge:
                if n not in nodeset:
                    nodeset.add(n)
                    pg.AddNode(n)

        for hyperedge in edges:
            if len(hyperedge) == 1: continue
            for i in range(0, len(hyperedge)-1):
                for j in range(i+1, len(hyperedge)):
                    i1, i2 = min(hyperedge[i], hyperedge[j]), max(hyperedge[i], hyperedge[j])
                    ret = pg.AddEdge(i1, i2)

        num_nodes = len(nodeset)
        effective_diameter = snap.GetBfsEffDiam(pg, num_nodes if num_nodes < 5000 else 5000, False)
        
        # Save Effective Diameter
        with open(outputdir + "effdiameter.txt", "w") as f:
            f.write(str(effective_diameter) + "\n")

def get_overlapness(inputpath, dataname, algorithmlist, portionlist, split_num, split_idx, recalculate_flag):
    print("Finding Overlapness...")

    if dataname is not None:
        _dataset = [dataname]
    else:
        _dataset = dataset

    if len(inputpath) > 0:
        dir_list_all = get_subdirectories([inputpath])
    else:
        dir_list_all = []
        # answer
        for dataname in _dataset:
            inputpath = "../../dataset/" + dataname + ".txt"
            outputdir = "../results/answer_dist/" + dataname + "/"
            dir_list_all.append((inputpath, outputdir))
        # algorithms
        for algoname in algorithmlist:
            for dataname in _dataset:
                for portion in portionlist:
                    # repeat
                    for i in range(1, repeat_time + 1):
                        dir_name = "../results/" + algoname + "/" + dataname + "_" + portion + "/" + str(i) + "/" 
                        if os.path.isdir(dir_name):
                            dir_list_all.append(dir_name)
                        # for mhname in mhlist:
                        #     mhdir_name = "../results/" + algoname + "_" + mhname + "/" + dataname + "_" + portion + "/" + str(i) + "/"
                        #     if os.path.isdir(mhdir_name):
                        #         dir_list_all.append(mhdir_name)
                    
                    dir_name = "../results/" + algoname + "/" + dataname + "_" + portion + "/"
                    if os.path.isdir(dir_name):
                        dir_list_all.append(dir_name)
                    # for mhname in mhlist:
                    #     mhdir_name = "../results/" + algoname + "_" + mhname + "/" + dataname + "_" + portion + "/"
                    #     if os.path.isdir(mhdir_name):
                    #         dir_list_all.append(mhdir_name)

        # split_gap = math.ceil(len(algorithmlist_nmh) / split_num)
        # end_idx = min(len(algorithmlist_nmh), split_gap * (split_idx + 1))
        # for algoname in algorithmlist_nmh[split_gap * split_idx : end_idx]:
        #     for dataname in _dataset:
        #         for portion in portionlist:
        #             for i in range(1, repeat_time + 1):
        #                 dir_name = "../results/" + algoname + "/" + dataname + "_" + portion + "/" + str(i) + "/" 
        #                 if os.path.isdir(dir_name):
        #                     dir_list_all.append(dir_name)
        #             dir_name = "../results/" + algoname + "/" + dataname + "_" + portion + "/"
        #             if os.path.isdir(dir_name):
        #                 dir_list_all.append(dir_name)

    split_gap = math.ceil(len(dir_list_all) / split_num)
    end_idx = min(len(dir_list_all), split_gap * (split_idx + 1))
    dir_list_all = dir_list_all[split_gap * (split_idx) : end_idx]            
    for dir_name in tqdm(dir_list_all):
        #print(dir_name)
        if type(dir_name) is tuple:
            inputpath, outputdir = dir_name
        else:
            inputpath = dir_name + "sampled_graph.txt"
            outputdir = dir_name

        flag_sampled_graph = os.path.isfile(inputpath)
        flag_ov = os.path.isfile(outputdir + "overlapness.txt")
        if flag_sampled_graph is False:
            continue
        elif (recalculate_flag is False) and (flag_ov is True):
            continue

        nodeset = set()
        sum_of_hyperedge_sizes  = 0

        with open(inputpath, "r") as f:
            for idx, line in enumerate(f.readlines()):
                line = line[:-1] # strip enter
                nodes = line.split(",")
                for i, node in enumerate(nodes):
                    if node not in nodeset:
                        nodeset.add(node)
                sum_of_hyperedge_sizes += len(nodes)
        overlapness = sum_of_hyperedge_sizes / len(nodeset)

        if sum_of_hyperedge_sizes == 0:
            print("Empty Input")
            return
        try:
            if os.path.isdir(outputdir) is False:
                continue
        except:
            print(inputpath)
            print(outputdir)
            print(dir_name)
            print(type(dir_name))

        if os.path.isdir(outputdir) is False:
                continue

        # Save overlapness
        with open(outputdir + "overlapness.txt", "w") as f:
            f.write(str(overlapness) + "\n")

def get_dist_from_dir(dir_path):
    if os.path.isfile(dir_path + "singular_values_full.txt"):
        s = []
        with open(dir_path + "singular_values_full.txt", "r") as f:
            for line in f.readlines():
                line = line[:-1]
                s.append(float(line))
        singular_values = {}
        singular_values["singular_value"] = s
    else:
        return -1
    
    if not os.path.isfile(dir_path + "sizewcc.txt"):
        return -1
    size_wcc = pd.read_csv(dir_path + "sizewcc.txt").sort_values(by=['size_wcc'])
    size_wcc['value'] = list(size_wcc['size_wcc'] / size_wcc['num_nodes'])
    size_wcc['size_wcc'] = list(range(1, len(size_wcc)+1))
    
    if not os.path.isfile(dir_path + "density.txt"):
        return -1
    density_pd = pd.read_csv(dir_path + "density.txt")
    density = float(density_pd["num edges"] / density_pd["num nodes"])

    if not os.path.isfile(dir_path + "overlapness.txt"):
        return -1
    with open(dir_path + "overlapness.txt", "r") as f:
        overlapness = float(f.readline())
    
    if not os.path.isfile(dir_path + "global_cc.txt"):
        return -1
    global_cc = -1
    with open(dir_path + "global_cc.txt", "r") as f:
        global_cc = f.readline()
        global_cc = float(global_cc)
    
    if not os.path.isfile(dir_path + "effdiameter.txt"):
        return -1
    effective_diam = -1
    with open(dir_path + "effdiameter.txt", "r") as f:
        effective_diam = f.readline()
        effective_diam = float(effective_diam)

    dist = {"singular_value" : singular_values, "size_wcc" : size_wcc,
            "density": density, "global_cc": global_cc, "effective_diameter": effective_diam, "overlapness" : overlapness}
    return dist

def save_all_evaluation(inputpath, dataname, algorithmlist, portionlist, split_num, split_idx):
    evallist = ["degree", "intersect", "pairdeg", "size", "Time", 
                "singular_value", "size_wcc", "global_cc", "density", "overlapness", "effective_diameter"]
                # "clustering_coef", "pathlength" , "nnz"

    answer_dir = "../results/answer_dist/"

    if dataname is not None:
        _dataset = [dataname]
    else:
        _dataset = dataset

    for dataname in _dataset:
        answer_dir = "../results/answer_dist/" + dataname + "/"
        answer_dist = get_dist_from_dir(answer_dir)

        if len(inputpath) > 0:
            dir_list_all = get_subdirectories([inputpath], dataname)
        else:
            dir_list_all = []
            for portion in portionlist:
                for algoname in algorithmlist:
                    for i in range(1, repeat_time + 1):
                        dir_name = "../results/" + algoname + "/" + dataname + "_" + portion + "/" + str(i) + "/" 
                        if os.path.isdir(dir_name):
                            dir_list_all.append(dir_name)
                    dir_name = "../results/" + algoname + "/" + dataname + "_" + portion + "/"
                    if os.path.isdir(dir_name):
                        dir_list_all.append(dir_name)
        split_gap = math.ceil(len(dir_list_all) / split_num)
        end_idx = min(len(dir_list_all), split_gap * (split_idx + 1))
        dir_list_all = dir_list_all[split_gap * (split_idx) : end_idx]               
        for dir_path in tqdm(dir_list_all, desc="( " + dataname + " )"):
            assert dataname in dir_path
            evaluation_result = {}
            if os.path.isfile(dir_path + "result.txt") is False:
                print(dir_path + " has no result.txt file")
                continue
            # Get dist
            data_dist = get_dist_from_dir(dir_path)
            if data_dist == -1:
                if os.path.isfile(dir_path + "entire_evaluation.txt"):
                    os.remove(dir_path + "entire_evaluation.txt")
                print("Error read result in ", dir_path)
                continue

            # When there exist all evaluations,
            with open(dir_path + "result.txt") as f:
                target_eval = ["degree", "intersect", "pairdeg", "size", "Time"]
                for line in f.readlines():
                    tmp = line[:-1].split(" : ")
                    if len(tmp) == 2:
                        ename, result = line[:-1].split(" : ")
                        if ename in target_eval:
                            evaluation_result[ename] = float(result)
                """
                if "grid_search" in dir_path:
                    agg_time = 0
                    tmp = dir_path.split("/")
                    #"../results/" + algoname + "/" + dataname + "_" + portion + "/" + str(i) + "/" 
                    tmp_algoname = tmp[2]
                    assert "grid_search" in tmp_algoname
                    tmp_algotype = tmp_algoname[:2]
                    prefix = "global_deg"
                    if "max" in dir_path:
                        prefix = "global_deg_max"
                    elif "avg" in dir_path:
                        prefix = "global_deg_avg"

                    agg_dir_list = [tmp_algotype + "/random_1.0000"]
                    agg_dir_list += [tmp_algotype + "/" + prefix + "_" + i for i in ["0.5000", "1.0000", "2.0000", "4.0000", "8.0000", "16.0000", "32.0000", "64.0000"]]
                    for algoname in agg_dir_list:
                        tmp[2] = algoname
                        new_dir_path = "/".join(tmp)
                        if os.path.isfile(new_dir_path + "result.txt"):
                            with open(new_dir_path + "result.txt") as f:
                                for line in f.readlines():
                                    ls = line[:-1].split(" : ")
                                    if len(ls) == 2:
                                        ename, result = line[:-1].split(" : ")
                                        if ename == "Time":
                                            agg_time += float(result)
                    evaluation_result["Time"] = agg_time
                """

            # Every time, write new entire_evaluation file!
            with open(dir_path + "entire_evaluation.txt", "w") as f:
                for evalname in evaluation_result.keys():
                    f.write(evalname + " : " + str(evaluation_result[evalname]) + "\n")

            assert "intersect" in evaluation_result, dir_path
            for evalname in evallist:
                if evalname not in evaluation_result:
                    if evalname in ["density", "overlapness", "global_cc", "effective_diameter"]:
                        # , "nnz"
                        # Difference
                        evaluation_result[evalname] = answer_dist[evalname] - data_dist[evalname]
                        evaluation_result[evalname + "_norm"] = (answer_dist[evalname] - data_dist[evalname]) / answer_dist[evalname]
                        with open(dir_path + "entire_evaluation.txt", "a") as f:
                            f.write(evalname + " : " + str(evaluation_result[evalname]) + "\n")
                            f.write(evalname + "_norm : " + str(evaluation_result[evalname + "_norm"]) + "\n")
                    elif evalname == "singular_value":
                        continue
                    elif answer_dist[evalname] is not None:
                        # D-Statistics
                        answer_dict = defaultdict(float)
                        data_dict = defaultdict(float)
                        keyset = set()
                        for idx in range(len(answer_dist[evalname][evalname])):
                            k = answer_dist[evalname][evalname][idx]
                            keyset.add(k)
                            v = answer_dist[evalname]["value"][idx]
                            answer_dict[k] = v
                        for idx in range(len(data_dist[evalname][evalname])):
                            k = data_dist[evalname][evalname][idx]
                            keyset.add(k)
                            v = data_dist[evalname]["value"][idx]
                            data_dict[k] = v
                        keys = [k for k in keyset]
                        keys = sorted(keys)
                        
                        answer_cumulsum = 0.0
                        data_cumulsum = 0.0
                        stat = 0.0
                        for k in keys:
                            if k in answer_dict:
                                assert answer_cumulsum <= answer_dict[k]
                                answer_cumulsum = answer_dict[k]
                                assert answer_cumulsum <= 1.0
                            if k in data_dict:
                                assert data_cumulsum <= data_dict[k]
                                data_cumulsum = data_dict[k]
                                assert data_cumulsum <= 1.0
                            if stat < math.fabs(answer_cumulsum - data_cumulsum):
                                stat = math.fabs(answer_cumulsum - data_cumulsum)
                        evaluation_result[evalname] = stat
                        with open(dir_path + "entire_evaluation.txt", "a") as f:
                            f.write(evalname + " : " + str(evaluation_result[evalname]) + "\n")   
            
            another_sv_eval_result_list = ["singular_value2"] #, "singular_value3"]
            for another_sv_result_path in another_sv_eval_result_list:
                if os.path.isfile(dir_path + another_sv_result_path + "_eval.txt"):
                    with open(dir_path + another_sv_result_path + "_eval.txt", "r") as rf:
                        stat2 = float(rf.readline())
                    with open(dir_path + "entire_evaluation.txt", "a") as f:
                        f.write("singular_value : " + str(stat2) + "\n")

def get_evaluation(dir_path):
    evallist = ["degree", "intersect", "pairdeg", "size", "Time",  
    "singular_value", "size_wcc", "global_cc", "density", "overlapness", "effective_diameter",
    "global_cc_norm", "density_norm", "overlapness_norm", "effective_diameter_norm"]
    # "nnz", "clustering_coef", "pathlength",
    
    evaluation = {}
    if os.path.isfile(dir_path + "/entire_evaluation.txt") is False:
        return -1
        
    with open(dir_path + "/entire_evaluation.txt" , "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.split("\n")[0]
            tmp = line.split(" : ")
            if len(tmp) == 2:
                eval_name, result = tmp
                if eval_name in evallist:
                    evaluation[eval_name] = float(result)

    for eval_name in evallist:
        if eval_name not in evaluation:
            return -1

    return evaluation

def aggregate_repeatition(inputpath, dataname, algorithmlist, portionlist, split_num, split_idx):
    evallist = ["degree", "intersect", "pairdeg", "size", "Time",  
    "singular_value", "size_wcc", "global_cc", "density", "overlapness", "effective_diameter",
    "global_cc_norm", "density_norm", "overlapness_norm", "effective_diameter_norm"]

    if dataname is not None:
        _dataset = [dataname]
    else:
        _dataset = dataset

    if len(inputpath) > 0:
        tmp = get_subdirectories([inputpath])
        dir_list_all = []
        for d in tmp:
            for i in range(1, repeat_time + 1):
                if d.endswith("/" + str(i) + "/"):
                    _d = d[:-2]
                    if _d not in dir_list_all:
                        dir_list_all.append(_d)
    else:
        dir_list_all = []
        for dataname in tqdm(_dataset):
            for portion in portionlist:
                for algoname in algorithmlist:
                    dir_list_all.append("../results/" + algoname + "/" + dataname + "_" + portion + "/")
    
    split_gap = math.ceil(len(dir_list_all) / split_num)
    end_idx = min(len(dir_list_all), split_gap * (split_idx + 1))
    dir_list_all = dir_list_all[split_gap * (split_idx) : end_idx]
    for d in tqdm(dir_list_all):
        eval_dict = defaultdict(list)
        arg_dict = {}
        min_dict = {}
        for i in range(1, repeat_time + 1):
            dirname = d + str(i)
            result = get_evaluation(dirname)
            if result == -1:
                print(dirname + "  FAIL")
                continue
            for ename in evallist:
                assert ename in result
            for ename in result.keys():
                if ename not in arg_dict:
                    arg_dict[ename] = i
                    min_dict[ename] = abs(result[ename])
                elif min_dict[ename] > abs(result[ename]):
                    arg_dict[ename] = i
                    min_dict[ename] = abs(result[ename])
                eval_dict[ename].append(result[ename])
        if len(eval_dict) == 0:
            print("No ", dirname)
            if os.path.isfile(d + "agg_entire_evaluation.txt"):
                os.remove(d + "agg_entire_evaluation.txt")
            if os.path.isfile(d + "agg_result_arg.txt"):
                os.remove(d + "agg_result_arg.txt")
            continue
        # print(eval_dict.keys())
        for ename in evallist:
            # if ename == "Time" and "search" in d:
            #     # read "time.txt":
            #     time_path = d + "time.txt"
            #     assert os.path.isfile(time_path)
            #     agg_time = 0
            #     with open(time_path) as f:
            #         agg_time = float(f.readline())
            #     assert agg_time != 0
            #     eval_dict["Time"] = agg_time
            #     arg_dict["Time"] = 0
            # else:                
            #     arg_dict[ename] = np.argmin( [abs(e) for e in eval_dict[ename]] ) + 1
            #     eval_dict[ename] = np.mean(eval_dict[ename])
            #arg_dict[ename] = np.argmin( [abs(e) for e in eval_dict[ename]] ) + 1
            eval_dict[ename] = np.mean(eval_dict[ename])


        with open(d + "agg_entire_evaluation.txt", "w") as r:
            for ename in eval_dict.keys():
                r.write(ename + " : " + str(eval_dict[ename]) + "\n")

        with open(d + "agg_result_arg.txt", "w") as r:
            for ename in arg_dict.keys():
                r.write(ename + " : " + str(arg_dict[ename]) + "\n")

if __name__ == "__main__":
    algorithmlist = ["es_grid_search", "es_skewness_search", 
                    "es_grid_search_ext", "es_skewness_search_ext",
                    "es_global_deg_avg_grid_search", "es_global_deg_max_grid_search",
                    "ns_grid_search", "ns/global_deg_1.0000", "es/global_deg_min_0.0000", "ns/global_deg_0.0000",
                    "tihs", "rw/rw_c_1", "ff/ff_c_0.51_0.20",
                    "mh/add_degree", "mh/exchange_degree", "mh/remove_degree",
                    "mh/add_avg", "mh/exchange_avg", "mh/remove_avg"]
                # "greedy/exchange/degree", "greedy/exchange/avg", 
    portionlist = ["0.30", "0.10", "0.20", "0.40", "0.50"]

    parser = argparse.ArgumentParser()
    parser.add_argument('--inputpath', required=False, default="")
    parser.add_argument('--read_inputs', required=False, action='store_true')
    parser.add_argument('--dataname', required=False)
    parser.add_argument('--algorithmlist', required=False, default="")
    parser.add_argument('--portionlist', required=False, default="")
    parser.add_argument('--portion', required=False, type=float)
    parser.add_argument('--diameter', required=False, action='store_true')
    parser.add_argument('--overlapness', required=False, action='store_true')
    parser.add_argument('--get_eval', required=False, action='store_true')
    parser.add_argument('--repeat_agg', required=False, action='store_true')
    parser.add_argument('--recalculate', required=False, action='store_true')
    parser.add_argument('--outputname', required=False, default='answer')
    parser.add_argument('--split_num', required=False, type=int, default=1)
    parser.add_argument('--split_idx', required=False, type=int, default=0)


    args = parser.parse_args()
    print("\n<DATA NAME>", args.dataname)

    if len(args.algorithmlist) > 0:
        algorithmlist = args.algorithmlist.split(",")
    if len(args.portionlist) > 0:
        portionlist = args.portionlist.split(",")

    if args.diameter:
        get_diameter(args.inputpath, args.dataname, algorithmlist, portionlist, args.split_num, args.split_idx, args.recalculate)
    if args.overlapness:
        get_overlapness(args.inputpath, args.dataname, algorithmlist, portionlist, args.split_num, args.split_idx, args.recalculate)

    if args.get_eval:
        save_all_evaluation(args.inputpath, args.dataname, algorithmlist, portionlist, args.split_num, args.split_idx)
    if args.repeat_agg:
        aggregate_repeatition(args.inputpath, args.dataname, algorithmlist, portionlist, args.split_num, args.split_idx)
