import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
from itertools import chain
import os
import math
from collections import defaultdict
import shutil
from copy import deepcopy
from scipy.stats import norm
import math
import matplotlib.pyplot as plt
import random
import argparse
plt.rcParams.update({'font.size': 20})

dataset = ["contact-high-school", "contact-primary-school", 
           "email-Enron-full", "email-Eu-full", 
           "tags-math-sx", "tags-ask-ubuntu" ,  
           "NDC-classes-full", "NDC-substances-full", 
           "coauth-MAG-Geology-full", "coauth-MAG-History-full", "threads-ask-ubuntu"]

portionlist = [0.1, 0.2, 0.3, 0.4, 0.5]

###################  Grid Search  #####################
def alpha_grid_search_all(algotype, opt, search_space):
    if "max" in opt:
        outputname = "maxdegree"
    elif "avg" in opt:
        outputname = "avgdegree"
    elif "min" in opt:
        outputname = "midas_grid_ablation"
    elif "ns" == algotype:
        outputname = "midas_ns"
        
    algorithmlist = []
    for a in search_space:
        if a == 0:
            if algotype == "es":
                algorithmlist.append(algotype + "/global_deg_min_0.0000")
            elif algotype == "ns":
                algorithmlist.append(algotype + "/global_deg_0.0000")
        else:
            algorithmlist.append(algotype + "/" + opt + "_%.4f" % (a))
    
    for portion in portionlist:
        portion_str = "%.2f" % (portion)
        for idx, data in enumerate(dataset):
            not_exist_flag = False
            best_performing_alpha = -1
            best_dstat = 1000
            for algoname in algorithmlist:
                dir_path = "../results/" + algoname + "/" + data + "_" + portion_str
                avg_deg = 0.0
                count = 0
                for i in range(1,4):
                    fname = dir_path + "/" + str(i) + "/result.txt"
                    if os.path.isfile(fname):
                        count += 1
                        with open(fname , "r") as f:
                            lines = f.readlines()
                            for line in lines:
                                line = line.split("\n")[0]
                                tmp = line.split(" : ")
                                if len(tmp) == 2:
                                    eval_name, result = tmp
                                    if eval_name == "degree":
                                        avg_deg += float(result)
                #assert count > 0, algoname + "/" + data
                if count == 0:
                    print(algoname, data, count)
                    not_exist_flag = True
                    break
                avg_deg /= count
                #result_ls.append(avg_deg)
                if avg_deg < best_dstat:
                    best_dstat = avg_deg
                    alpha_str = algoname.split("_")[-1]
                    best_performing_alpha =  float(alpha_str)

                if count < 3:
                    print(algoname, data, count)
                    
            if not_exist_flag is True:
                best_performing_alpha = -1
            else:
                # Save Result
                if os.path.isdir("../results/" + outputname + "/") is False:
                    os.makedirs("../results/" + outputname + "/")
                with open("../results/" + outputname + "/search_result.txt", "a+") as f:
                    if algotype == "ns":
                        path = "../results/" + algotype + "/global_deg_%.4f" % (best_performing_alpha)
                    elif best_performing_alpha == 0:
                        path = "../results/" + algotype + "/global_deg_min_%.4f" % (best_performing_alpha)
                    else:
                        path = "../results/" + algotype + "/" + opt + "_%.4f" % (best_performing_alpha)
                    f.write(path + "/" + data + "_" + portion_str + "\n")
                        
################### Observation 1 #####################
def get_dirpath(dir_path):
    _dir_path = dir_path
    
    if "answer" in dir_path:
        return _dir_path
    
    pick = random.randrange(1,4)
    _dir_path = dir_path + "/" + str(pick)
                    
    return _dir_path

def get_dist(dir_path):
    dir_path = get_dirpath(dir_path)
        
    if not os.path.isfile(dir_path + "/" + "degree_dist.txt"):
        print("No ", dir_path + "/" + "degree_dist.txt")
        return -1
        
    deg = pd.read_csv(dir_path +  "/degree_dist.txt").sort_values(by=['degree'])
    deg['value'] = deg['value'].cumsum()
    
    dist = {'degree': deg}
    return dist

def plot_degree_figure(algorithm_dist,dataname, portion, color_dict, line_dict, dirname):
    plt.figure(figsize=(5.5,3.4), dpi=120)
    eval_name = "degree"
    
    for algo in algorithm_dist.keys():
        color, line = color_dict[algo], line_dict[algo]
        if algo == "answer":
            plt.plot(algorithm_dist[algo][eval_name][eval_name], algorithm_dist[algo][eval_name]['value'], color=color, alpha=1.0, linewidth=12, linestyle=line, label=algo)
        else:
            plt.plot(algorithm_dist[algo][eval_name][eval_name], algorithm_dist[algo][eval_name]['value'], color=color, alpha=1.0, linewidth=4, linestyle=line, label=algo)
        
    
    plt.xscale('log')
    plt.xlabel("Degree", fontsize=20)
    plt.ylabel("Cumulative\nProbability", fontsize=20)
    
    ax = plt.gca()
    ax.tick_params(labelcolor='#4B4B4B', labelsize=16)
    yticks = ax.get_yticks()[1:-1]
    plt.tight_layout()
    # plt.legend()
    
    savedir = "figures/Appendix/" + dirname + "/%.1f/"% (portion)
    if os.path.isdir(savedir) is False: 
        os.makedirs(savedir)
    savename = savedir + dataname + ".jpg"
    plt.savefig(savename, bbox_inches='tight')
    plt.show()
    plt.close()

def observation(algotype, opt, search_space, dataname, portion):
    if "max" in opt:
        outputname = "maxdegree"
    elif "avg" in opt:
        outputname = "avgdegree"
    elif "min" in opt:
        outputname = "midas_grid_ablation"
    elif "ns" == algotype:
        outputname = "midas_ns"
        
    algorithmlist = []
    for a in search_space:
        if a == 0:
            if algotype == "es":
                algorithmlist.append(algotype + "/global_deg_min_0.0000")
            elif algotype == "ns":
                algorithmlist.append(algotype + "/global_deg_0.0000")
        else:
            algorithmlist.append(algotype + "/" + opt + "_%.4f" % (a))

    colors = ["#FFC81E", "#80E12A", "#64CD3C", "#228B22", "#147814", "#0A8A8A", "#00BFFF", "#1E90FF", "#0064FF", "#8572EE"]
    color_dict ={"answer" : "black"}
    line_dict = {"answer" : "dashdot"}
    for idx, algo in enumerate(algorithmlist):
        color_dict[algo] = colors[idx]
        line_dict[algo] = "solid"

    algo_dist = {}
    algo_dist["answer"] = get_dist("../results/answer_dist/" + dataname)
    for algo_name in algorithmlist:
        algo_dir = "../results/" + algo_name + "/" + dataname + "_%.2f" % (portion)
        ret = get_dist(algo_dir)
        if ret == -1:
            continue
        algo_dist[algo_name] = ret
    plot_degree_figure(algo_dist, dataname, portion, color_dict, line_dict, "obs1/" + outputname )

################### Degree Distribution #####################
def compare_degree_dist(dataname, portion):
    algorithmlist=["midas_grid_ablation", "maxdegree", "avgdegree", "midas_ns"]
    algo_color = {
        "midas_grid_ablation": "#377eb8",
        "maxdegree": "#4daf4a",
        "avgdegree": "#984ea3",
        "midas_ns": "#e41a1c",
        "answer" : "black"
    }
    line_dict = {"maxdegree" : "solid", "avgdegree" : "solid", "midas_ns" : "solid", "midas_grid_ablation" : "solid" , "answer" : "dashdot"}
    
    algo_dist = {}
    algo_dist["answer"] = get_dist("../results/answer_dist/" + dataname)
    for algo_name in algorithmlist:
        algo_dir = "../results/" + algo_name + "/" + dataname + "_%.2f" % (portion)
        ret = get_dist(algo_dir)
        if ret == -1:
            print(algo_dir)
            continue
        algo_dist[algo_name] = ret
    plot_degree_figure(algo_dist, dataname, portion, algo_color, line_dict, "degree_dist")
    

################### According to Portion #####################
def analyze_evaluation():
    algorithmlist=["midas_grid_ablation", "maxdegree", "avgdegree", "midas_ns"]
    color_dict = {
        "midas_grid_ablation": "#377eb8",
        "maxdegree": "#4daf4a",
        "avgdegree": "#984ea3",
        "midas_ns": "#e41a1c",
        "answer" : "black"
    }
    marker_dict = {
        "midas_grid_ablation": "o",
        "maxdegree": "o",
        "avgdegree": "o",
        "midas_ns": "o",
        "answer" : ","
    }
    algo2portion_ranklist = defaultdict(list)
    algo2portion_normlist = defaultdict(list)
    algo2portion_degreedstat = defaultdict(list)
    
    for portion in portionlist:
        portion_str = "%.2f" % (portion)
        d = pd.read_csv("./csvs/" + portion_str + "_ablation_agg.csv")
        for i, row in d.iterrows():
            algoname = row["algorithm"]
            if row["algo opt"] != "-":
                algoname += "/" + row["algo opt"]
            if row["eval opt"] != "-":
                algoname += "/" + row["eval opt"]
            algo2portion_ranklist[algoname].append(float(row["aggregate rank"]))
            algo2portion_degreedstat[algoname].append(float(row["degree avg"]))

        d = pd.read_csv("./csvs/" + portion_str + "_ablation_norm_agg.csv")
        d["algoname"] = d["algorithm"] + "/" + d["algo opt"]
        for i, row in d.iterrows():
            algoname = row["algorithm"]
            if row["algo opt"] != "-":
                algoname += "/" + row["algo opt"]
            if row["eval opt"] != "-":
                algoname += "/" + row["eval opt"]
            algo2portion_normlist[algoname].append(float(row["avg"]))

    # Plot Ranking vs. Portion
    plt.figure(figsize=(4.2, 4.0), dpi=120)
    portionlabel = [10, 20, 30, 40, 50]
    for algo in algorithmlist:
        if "midas_grid_ablation" in algo:
            plt.plot(portionlabel, algo2portion_ranklist[algo], marker=marker_dict[algo], markersize=10, label=algo, c=color_dict[algo],alpha=1.0, linewidth=8)
        else:
            plt.plot(portionlabel, algo2portion_ranklist[algo], marker=marker_dict[algo], markersize=10, label=algo, c=color_dict[algo], alpha=1.0, linewidth=8)

    plt.ylabel("Avg. Rank", fontsize=20)
    plt.xlabel("Sampling Portion (%)", fontsize=20)
    ax = plt.gca()
    ax.set_xticks(portionlabel)
    ax.tick_params(labelcolor='#4B4B4B', labelsize=18)
    # plt.legend(bbox_to_anchor=(-0.2, 1), fontsize=15)
    plt.tight_layout()
    savedir = "figures/Appendix/"
    if os.path.isdir(savedir) is False: 
        os.makedirs(savedir)
    savename = savedir + "ranking.jpg"
    plt.savefig(savename, bbox_inches='tight')
    plt.close()
    
    # Plot Z-Score vs. Portion
    plt.figure(figsize=(4.45, 4), dpi=120)
    for algo in algorithmlist:
        if "midas_grid_ablation" in algo:
            plt.plot(portionlabel, algo2portion_normlist[algo], marker=marker_dict[algo],  markersize=10, label=algo, c=color_dict[algo],alpha=1.0, linewidth=8)
        else:
            plt.plot(portionlabel, algo2portion_normlist[algo], marker=marker_dict[algo],   markersize=10, label=algo, c=color_dict[algo], alpha=1.0, linewidth=8)

    plt.ylabel("Avg. Z-Score", fontsize=20)
    plt.xlabel("Sampling Portion (%)", fontsize=20)
    ax = plt.gca()
    ax.set_xticks(portionlabel)
    ax.tick_params(labelcolor='#4B4B4B', labelsize=18)
    # plt.legend(bbox_to_anchor=(-0.2, 1), fontsize=15)
    ax.yaxis.set_label_coords(-0.25,0.48)

    plt.tight_layout()
    savedir = "figures/Appendix/"
    if os.path.isdir(savedir) is False: 
        os.makedirs(savedir)
    savename = savedir + "zscore.jpg"
    plt.savefig(savename, bbox_inches='tight')
    plt.close()
    
    # Plot D-Statistics in Degree vs. Portion
    plt.figure(figsize=(4.6, 4), dpi=120)
    portionlabel = [10,20,30,40,50]
    for algo in algorithmlist:
        if "midas_grid_ablation" in algo:
            plt.plot(portionlabel, algo2portion_degreedstat[algo], label=algo, alpha=1.0, color=color_dict[algo], linewidth=8, marker="o",  markersize=10)
        else:
            plt.plot(portionlabel, algo2portion_degreedstat[algo], label=algo, alpha=1.0, color=color_dict[algo], linewidth=8, marker="o",  markersize=10)
    ax = plt.gca()
    ax.tick_params(labelcolor='#4B4B4B', labelsize=18)
    ax.set_xticks(portionlabel)
    # ax.tick_params(axis='x', labelrotation=90, labelsize=20)

    plt.ylabel("D-Statistics in\nDegree Distribution", fontsize=20)
    plt.xlabel("Sampling Portion (%)", fontsize=20)

    # plt.title( r"$\bf{Degree Dstat.}$", fontsize=18)
    plt.tight_layout()
    savedir = "figures/Appendix/"
    if os.path.isdir(savedir) is False: 
        os.makedirs(savedir)
    savename = savedir + "dstat.jpg"
    plt.savefig(savename, bbox_inches='tight')
    print(savename)
    
if __name__ == "__main__":
    search_space = [0, 0.5 ,1, 2, 4, 8, 16, 32, 64]
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--select', required=True, type=int, help="0: Alpha Grid Search, 1: Analyze on Observation1, 2: Degree Distribution, 3: Evaluation after run 'analyze_result.py' ")
    parser.add_argument('--algotype', required=False, type=str, help="(ns/es) node selection or hyperedge selection")
    parser.add_argument('--opt', required=False, type=str, help="(global_deg_min/global_deg_max/global_deg_avg/global_deg) for midas_grid_ablation / max degree / avg degree / midas_ns")
    parser.add_argument('--data', required=False, default="email-Eu-full", type=str)
    parser.add_argument('--portion', required=False, default=0.3, type=float)
    args = parser.parse_args()
    
    if args.select == 0:
        alpha_grid_search_all(args.algotype, args.opt, search_space)
    elif args.select == 1:
        observation(args.algotype, args.opt, search_space, args.data, args.portion)
    elif args.select == 2:
        compare_degree_dist(args.data, args.portion)
    elif args.select == 3:
        analyze_evaluation()
        