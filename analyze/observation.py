import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import defaultdict
import seaborn as sns
import math
import argparse
import os
from sklearn.linear_model import LinearRegression

dataset = ["email-Enron-full", "email-Eu-full",
           "contact-high-school", "contact-primary-school",
           "NDC-classes-full", "NDC-substances-full",
           "tags-ask-ubuntu", "tags-math-sx", "threads-ask-ubuntu",
           "coauth-MAG-History-full", "coauth-MAG-Geology-full"]

portionlist = ["0.30", "0.20", "0.10", "0.40", "0.50"]

###################### Function ######################
def correlation_matrix():
    #  --------------   Rank   --------------
    target_cols = ['rank degree',
                    'rank size',
                    'rank pairdeg',
                    'rank intersect',
                    'rank singular_value',
                    'rank size_wcc',
                    'rank global_cc_norm',
                    'rank density_norm',
                    'rank overlapness_norm',
                    'rank effective_diameter_norm',
                    'eval avg rank'
                    ]

    label_col = {'rank degree': 'Degree',
                'rank intersect': 'Int. Size',
                'rank pairdeg': 'Pair Degree',
                'rank size': 'Size',
                'rank size_wcc': 'CC',
                'rank global_cc_norm': 'GCC',
                'rank density_norm': 'Density',
                'rank overlapness_norm': 'Overlapness',
                'rank effective_diameter_norm': 'Diameter',
                'rank singular_value': 'SV',
                'eval avg rank': 'Avg. Ranking'
    }
    all_dict = defaultdict(list)
    for portion_str in portionlist:
        for didx, dataname in enumerate(dataset):
            d = pd.read_csv("./csvs/" + dataname + "/" + portion_str + "_baseline.csv")
            for col in target_cols:
                all_dict[label_col[col]] = all_dict[label_col[col]] + list(d[col])
    d = pd.DataFrame.from_dict(all_dict)
    corrM = d.corr()
    fig, ax = plt.subplots( figsize=(13,12) )
    mask = np.zeros_like(corrM, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True
    sns.heatmap(corrM, 
            cmap = 'RdYlBu_r', cbar=False,
            vmin = -1,vmax = 1
           )
    names = list(d.columns)
    ymin, ymax = ax.get_ylim()
    ax.text(-3.6, ymin + 0.7, '/ Z-Score', color="black", fontsize = 50)

    plt.yticks(np.arange(0.5, len(names), 1), names, fontsize=50, rotation=0)
    ax.set(xticklabels=[])
    ax.set(xlabel=None)
    ax.tick_params(bottom=False)
    savedir = "figures/Correlation/"
    if os.path.isdir(savedir) is False:
        os.makedirs(savedir)
    savefname = savedir + 'correlation_baseline_ranking.jpg'
    plt.savefig(savefname, bbox_inches='tight')
    plt.close()

    #  --------------   Z-Score   --------------
    target_cols = ['degree',
                    'size',
                    'pairdeg',
                    'intersect',
                    'singular_value',
                    'size_wcc',
                    'global_cc_norm',
                    'density_norm',
                    'overlapness_norm',
                    'effective_diameter_norm',
                    'avg'
    ]

    label_col = {'degree': 'Degree',
                'intersect': 'Int. Size',
                'pairdeg': 'Pair Degree',
                'size': 'Size',
                'size_wcc': 'CC',
                'global_cc_norm': 'GCC',
                'density_norm': 'Density',
                'overlapness_norm': 'Overlapness',
                'effective_diameter_norm': 'Diameter',
                'singular_value': 'SV',
                'avg': 'AVG. Z-Score'
    }
    all_dict = defaultdict(list)
    for portion_str in portionlist:
        for didx, dataname in enumerate(dataset):
            d = pd.read_csv("./csvs/" + dataname + "/" + portion_str + "_baseline_norm.csv")
            for col in target_cols:
                all_dict[label_col[col]] = all_dict[label_col[col]] + list(d[col])

    d = pd.DataFrame.from_dict(all_dict)
    corrM = d.corr()

    fig, ax = plt.subplots( figsize=(16,12) )
    mask = np.zeros_like(corrM, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True

    sns.heatmap(corrM, 
                cmap = 'RdYlBu_r',
                cbar_kws={"ticks": [-1.0, -0.5, 0.0, 0.5, 1.0]},
                vmin = -1,vmax = 1
            )

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=58)
    names = list(d.columns)
    plt.yticks(np.arange(0.5, len(names), 1), names, fontsize=50)
    ax.set(xticklabels=[])
    ax.set(yticklabels=[])
    ax.set(xlabel=None)
    ax.set(ylabel=None)
    ax.tick_params(left=False)
    ax.tick_params(bottom=False)
    savedir = "figures/Correlation/"
    if os.path.isdir(savedir) is False:
        os.makedirs(savedir)
    savefname = savedir + 'correlation_baseline_z.jpg'
    plt.savefig(savefname, bbox_inches='tight')
    plt.close()

###################### Function ######################
def get_dirpath(dir_path):
    if "answer" in dir_path:
        return dir_path
    
    _dir_path = dir_path
    # assert os.path.isfile(dir_path + "/agg_result_arg.txt"), dir_path
    if os.path.isfile(dir_path + "/agg_result_arg.txt") is False:
        _dir_path = dir_path + "/1/"
        print(dir_path)
    else: 
        with open(dir_path + "/agg_result_arg.txt", "r") as f:
            for line in f.readlines():
                ename, arg = line[:-1].split(" : ")
                if ename == "degree":
                    _dir_path = dir_path + "/" + arg + "/"

    return _dir_path

def get_avgdeg_density(_dir_path):
    dir_path = get_dirpath(_dir_path)
    
    deg = {}
    _deg = pd.read_csv(dir_path +  "degree_dist.txt").sort_values(by=['degree'])
    deg['value'] = list(_deg['value']) # not cumsum
    deg['degree'] = list(_deg['degree'])
    
    avg_deg = 0
    for i in range(len(deg["degree"])):
        avg_deg += (deg['degree'][i] * deg['value'][i])
        assert deg['value'][i] <= 1, str(i) + "/" + str(len(deg['value'])) + "/" + str(deg['value'][i])
        
    density = -1
    with open(dir_path + "density.txt", "r") as f:
        f.readline()
        _V, _E = f.readline().split(",")
        V, E = int(_V), int(_E)
        density = E / V
    
    return avg_deg, density

def get_avgdeg_overlapness(_dir_path):
    dir_path = get_dirpath(_dir_path)
    
    deg = {}
    _deg = pd.read_csv(dir_path +  "degree_dist.txt").sort_values(by=['degree'])
    deg['value'] = list(_deg['value']) # not cumsum
    deg['degree'] = list(_deg['degree'])
    
    avg_deg = 0
    for i in range(len(deg["degree"])):
        avg_deg += (deg['degree'][i] * deg['value'][i])
        assert deg['value'][i] <= 1, str(i) + "/" + str(len(deg['value'])) + "/" + str(deg['value'][i])
        
    overlapness = -1
    with open(dir_path + "overlapness.txt", "r") as f:
        overlapness = f.readline()
        overlapness = float(overlapness)
    
    return avg_deg, overlapness

def correlation_with_degree(portion):
    portion_str = "%.2f" % (portion)
    ################### STYLE ######################
    color_dict = {
        # NS
        "ns/global_deg_0.0000": "#e41a1c",
        "ns/global_deg_1.0000": "#FF7A85", #"#F56E6E",
        "rw/rw_c_1": "#ff7f00",
        "ff/ff_c_0.51_0.20": "#FFC81E", #"#984ea3",
        
        # ES
        "es/global_deg_min_0.0000": "#377eb8",
        "tihs": "#4daf4a",
        
        # ANSWER    
        "norm" : "grey",
        'answer': "#000000" 
    }
    marker_dict = {
        # NS
        "ns/global_deg_0.0000": ">", # red
        "ns/global_deg_1.0000": "^", # red
        # ES
        "es/global_deg_min_0.0000": "D", # green
        "tihs": "P", # green
        # Topological
        "rw/rw_c_1": "s",
        "ff/ff_c_0.51_0.20": "<",
        
        # ANSWER    
        "norm" : ",",
        'answer': "," 
    }
    data_labeling_dic = {
        "email-Enron-full" : "Email\ Enron",
        "email-Eu-full" : "Email\ Eu",
        "NDC-classes-full" : "NDC\ Class",
        "NDC-substances-full" : "NDC\ Substance",
        "contact-high-school" : "Contact\ High",
        "contact-primary-school" : "Contact\ Primary",
        "tags-ask-ubuntu" : "Tags\ Ubuntu",
        "tags-math-sx" : "Tags\ Math",
        "threads-ask-ubuntu" : "Threads\ Ubuntu",
        "threads-math-sx" : "Threads\ Math",
        "coauth-MAG-Geology-full" : "Coauth\ Geology",
        "coauth-MAG-History-full" : "Coauth\ History",
        "coauth-DBLP-full" : "Coauth\ DBLP",
    }
    #################################################

    algorithmlist=["es/global_deg_min_0.0000", "tihs", 
               "ns/global_deg_0.0000", "ns/global_deg_1.0000", 
               "rw/rw_c_1", "ff/ff_c_0.51_0.20"]

    #  ----------------   Density   ----------------
    for idx, dataname in enumerate(dataset):
        data_label = data_labeling_dic[dataname]

        X, Y = [], []
        ans_avgdeg, ans_density = get_avgdeg_density("../results/answer_dist/" + dataname + "/")
        
        avgdeg_dict = {}
        for algo_name in algorithmlist:
            dir_path = "../results/" + algo_name
            algo_dir = dir_path + "/" + dataname + "_" + portion_str + "/"
            avgdeg, density = get_avgdeg_density(algo_dir)
            Y.append(density)
            X.append(avgdeg)
        
        # Plotting
        plt.figure(figsize=(3.5,2.6), dpi=120)
        xmin, xmax = min(min(X)*0.9, ans_avgdeg), max(max(X)*1.1, ans_avgdeg)
        ymin, ymax = min(min(Y)*0.8, ans_density), max(max(Y)*1.2, ans_density)
        plt.scatter(ans_avgdeg, ans_density, c = "black", s = 500, alpha=0.9)
        xs, ys = [], []
        for aidx, algo in enumerate(algorithmlist):
            plt.scatter(X[aidx], Y[aidx], marker=marker_dict[algo], c = color_dict[algo], s = 500, alpha=0.9, zorder=2)
            xs.append(X[aidx])
            ys.append(Y[aidx])
        X, Y = np.array(xs).reshape(-1, 1), np.array(ys).reshape(-1, 1)
        reg = LinearRegression().fit(X, Y)
        coef, intercept = reg.coef_, reg.intercept_
        xmin, xmax = min(min(xs), ans_avgdeg), max(max(xs), ans_avgdeg)
        X = np.arange(xmin * 0.7, xmax * 1.2).reshape(-1,1)
        plt.plot(X, reg.predict(X), alpha=0.5, color="black", linewidth=3, linestyle="dashed", zorder=1)
        plt.xlabel('Avg. degree', fontsize=20)
        plt.ylabel('Density', fontsize=20)
        ax = plt.gca()
        ax.tick_params(labelcolor='#4B4B4B', labelsize=16)
        cor = np.corrcoef(np.array(xs), np.array(ys))[0,1]
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        plt.ylim(ymin*0.7, ymax*1.1)
        ax.text(0.4, 0.95, "corr.=%.2f" % (cor), transform=ax.transAxes, fontsize=18, verticalalignment='top')
        plt.tight_layout()
        savedir = "figures/Correlation/degree_density/%.1f" % (float(portion_str)) +"/"
        if os.path.isdir(savedir) is False:
            os.makedirs(savedir)
        savename = savedir + dataname + ".jpg"
        plt.savefig(savename, bbox_inches='tight')
        plt.show()
        plt.close()
        
    #  --------------   Overlapness   --------------
    for idx, dataname in enumerate(dataset):
        data_label = data_labeling_dic[dataname]
        X, Y = [], []
        ans_avgdeg, ans_overlapness = get_avgdeg_overlapness("../results/answer_dist/" + dataname + "/")
        for algo_name in algorithmlist:
            dir_path = "../results/" + algo_name
            algo_dir = dir_path + "/" + dataname + "_" + portion_str + "/"
            avgdeg, overlapness = get_avgdeg_overlapness(algo_dir)
            Y.append(overlapness)
            X.append(avgdeg)

        plt.figure(figsize=(3.3,2.6), dpi=120)
        xmin, xmax = min(min(X)*0.9, ans_avgdeg), max(max(X)*1.1, ans_avgdeg)
        ymin, ymax = min(min(Y)*0.8, ans_overlapness), max(max(Y)*1.2, ans_overlapness)
        plt.scatter(ans_avgdeg, ans_overlapness, color="black", alpha=0.9, s=500)
        xs, ys = [], []
        for aidx, algo in enumerate(algorithmlist):
            plt.scatter(X[aidx], Y[aidx], marker=marker_dict[algo], c = color_dict[algo], s = 500, alpha=1.0, zorder=2)
            xs.append(X[aidx])
            ys.append(Y[aidx])
  
        X, Y = np.array(xs).reshape(-1, 1), np.array(ys).reshape(-1, 1)
        reg = LinearRegression().fit(X, Y)
        coef, intercept = reg.coef_, reg.intercept_
        xmin, xmax = min(min(xs), ans_avgdeg), max(max(xs), ans_avgdeg)
        X = np.arange(xmin * 0.7, xmax * 1.2).reshape(-1,1)
        plt.plot(X, reg.predict(X), alpha=0.5, color="black", linewidth=3, linestyle="dashed", zorder=1)
        plt.xlabel('Avg. degree', fontsize=20)
        plt.ylabel('Overlapness', fontsize=20)
        ax = plt.gca()
        ax.tick_params(labelcolor='#4B4B4B', labelsize=16)
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        cor = np.corrcoef(np.array(xs), np.array(ys))[0,1]
        ax.text(0.4, 0.95, "corr.=%.2f" % (cor), transform=ax.transAxes, fontsize=18, verticalalignment='top')
        plt.tight_layout()
        savedir = "figures/Correlation/degree_overlap/%.1f" % (float(portion_str)) +"/"
        if os.path.isdir(savedir) is False:
            os.makedirs(savedir)
        savename = savedir + dataname + ".jpg"
        plt.savefig(savename, bbox_inches='tight')
        plt.show()
        plt.close()

###################### Function ######################
def get_dist(dir_path, portion=None):
    dir_path = get_dirpath(dir_path)
        
    for name in ['intersect', 'pairdeg', 'size', 'degree']:
        if not os.path.isfile(dir_path + "/" + name + "_dist.txt"):
            print("No ", dir_path + "/" + name + "_dist.txt")
            return -1
        
    its = pd.read_csv(dir_path + "/intersect_dist.txt").sort_values(by=['intersect'])
    its['value'] = its['value'].cumsum()
    pairdeg = pd.read_csv(dir_path + "/pairdeg_dist.txt").sort_values(by=['pairdeg'])
    pairdeg['value'] = pairdeg['value'].cumsum()
    size = pd.read_csv(dir_path + "/size_dist.txt").sort_values(by=['size'])
    size['value'] = size['value'].cumsum()
    deg = pd.read_csv(dir_path +  "/degree_dist.txt").sort_values(by=['degree'])
    deg['value'] = deg['value'].cumsum()
    
    dist = {'intersect': its, 'pairdeg': pairdeg, 'size': size, 'degree': deg}
    
    if portion is None: #original
        avg_helper = []
        with open(dir_path + "/result.txt" , "r") as f:
            lines = f.readlines()
            for line in lines:
                line = line.split("\n")[0]
                tmp = line.split(" : ")
                if len(tmp) == 2:
                    eval_name, result = tmp
                    if eval_name in ["degree", "intersect", "pairdeg", "size"]:
                        avg_helper.append(float(result))
                        dist[eval_name + "_dstat"] = float(result)
                    elif eval_name == "Time":
                        dist[eval_name] = float(result)
        dist["avg"] = np.mean(avg_helper)
    else:
        for eval_name in ["degree", "intersect", "pairdeg", "size"]:
            dist[eval_name + "_dstat"] = 0.0
        dist['Time'] = 0.0
        dist["avg"] = 0.0
    
    return dist

def midas_observation1(dataname, portion):
    algorithmlist = ["es/global_deg_min_0.0000", "es/global_deg_min_0.5000", "es/global_deg_min_1.0000", "es/global_deg_min_2.0000", "es/global_deg_min_4.0000", "es/global_deg_min_8.0000",  "es/global_deg_min_16.0000", "es/global_deg_min_32.0000","es/global_deg_min_64.0000"]

    ################### STYLE ######################
    colors = ["#FFC81E", "#80E12A", "#64CD3C", "#228B22", "#147814", "#0A8A8A", "#00BFFF", "#1E90FF", "#0064FF", "#8572EE"]
    color_dict = { "es/global_deg_min_0.0000": "#FFC81E",
                    "es/global_deg_min_0.5000": "#80E12A",
                    "es/global_deg_min_1.0000": "#64CD3C",
                    "es/global_deg_min_2.0000": "#228B22",
                    "es/global_deg_min_4.0000": "#147814",  
                    "es/global_deg_min_8.0000": "#0A8A8A",
                    "es/global_deg_min_16.0000": "#00BFFF",
                    "es/global_deg_min_32.0000": "#0078FF", #"#1E90FF",
                    "es/global_deg_min_64.0000": "#5A5AFF",
                    "answer": "black"
    }
    line_dict = {   "es/global_deg_min_0.0000": "solid",
                    "es/global_deg_min_0.5000": "solid",
                    "es/global_deg_min_1.0000": "solid",
                    "es/global_deg_min_2.0000": "solid",
                    "es/global_deg_min_4.0000": "solid",  
                    "es/global_deg_min_8.0000": "solid",
                    "es/global_deg_min_16.0000": "solid",
                    "es/global_deg_min_32.0000": "solid",
                    "es/global_deg_min_64.0000": "solid",
                    "answer": "dashdot"
    }
    marker_dict = { "es/global_deg_min_0.0000": ",",
                    "es/global_deg_min_0.5000": ",",
                    "es/global_deg_min_1.0000": ",",
                    "es/global_deg_min_2.0000": ",",
                    "es/global_deg_min_4.0000": ",",  
                    "es/global_deg_min_8.0000": ",",
                    "es/global_deg_min_16.0000": ",",
                    "es/global_deg_min_32.0000": ",",
                    "es/global_deg_min_64.0000": ",",
                    "answer": ","
    }
    ################################################
    # Construct algorithm_dist
    algorithm_dist = {}
    algorithm_dist["answer"] = get_dist("../results/answer_dist/" + dataname, portion)
    for algo_name in algorithmlist:
        algo_dir = "../results/" + algo_name + "/" + dataname + "_%.2f" % (portion)
        ret = get_dist(algo_dir)
        if ret == -1:
            continue
        algorithm_dist[algo_name] = ret
    
    plt.figure(figsize=(5.5,3.3), dpi=120)
    eval_name = "degree"
    for algoname in algorithm_dist.keys():
        color, line, marker= color_dict[algoname], line_dict[algoname], marker_dict[algoname]
        if algoname == "norm":
            plt.plot(algorithm_dist[algoname][eval_name][eval_name], algorithm_dist[algoname][eval_name]['value'], color=color, alpha=0.3, linewidth=15, linestyle=line, label=algoname, marker=marker, markersize=10)
        elif algoname == "answer":
            plt.plot(algorithm_dist[algoname][eval_name][eval_name], algorithm_dist[algoname][eval_name]['value'], color=color, alpha=1.0, linewidth=18, linestyle=(0,(5,1)), label=algoname, zorder=1) #, marker=marker, markersize=10)
        else:
            plt.plot(algorithm_dist[algoname][eval_name][eval_name], algorithm_dist[algoname][eval_name]['value'], color=color, alpha=1.0, linewidth=4, linestyle="solid", label=algoname) #, marker=marker, markersize=5)
    plt.xscale('log')
    plt.xlabel("Degree", fontsize=20)
    plt.ylabel("Cumulative\nProbability", fontsize=20)
    
    ax = plt.gca()
    ax.tick_params(labelcolor='#4B4B4B', labelsize=16)
    yticks = ax.get_yticks()[1:-1]
    left_flag = False
    for y in yticks:
        if len(str(y)) > 3:
            left_flag = True
            break
    locmaj = mpl.ticker.LogLocator(numticks=5)
    ax.xaxis.set_major_locator(locmaj)

    if left_flag is False:
        ax.yaxis.set_label_coords(-0.3,0.48)
    else:
        ax.yaxis.set_label_coords(-0.3,0.48)
    
    plt.tight_layout()
    savedir = "figures/Observation/%.1f"% (portion) + "/dist/"
    if os.path.isdir(savedir) is False:
        os.makedirs(savedir)
    savefname = savedir + dataname + ".jpg"
    plt.savefig(savefname, bbox_inches='tight')
    plt.show()


###################### Function ######################
def get_avg_degreeDstat(algoname, data, portion_str):
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
    assert count > 0, algoname + "/" + data + "/" + portion_str
    avg_deg /= count

    return avg_deg

def midas_observation2():
    algorithmlist = ["es/global_deg_min_%.4f" % (2 ** i) for i in np.arange(-3,6.5,0.5)]
    algorithmlist += ["es/global_deg_min_0.0000"]
    skewness = {
        "email-Enron-full": 1.220866,
        "email-Eu-full" : 2.521832,
        "contact-high-school" : 0.480156,
        "contact-primary-school": 0.314277,
        "NDC-classes-full" : 8.603801,
        "NDC-substances-full" : 8.702568,
        "tags-ask-ubuntu" : 10.307134,
        "tags-math-sx" : 6.803574,
        "threads-ask-ubuntu" : 52.012624,
        "threads-math-sx" : 57.982697,
        'coauth-DBLP-full' : 14.744422,
        "coauth-MAG-History-full" : 52.414521,
        "coauth-MAG-Geology-full" : 12.639768,
    }
    dataset_color = {
        "email-Enron-full": "#e41a1c",
        "email-Eu-full": "#e41a1c",
        "contact-high-school": "#377eb8",
        "contact-primary-school": "#377eb8",
        "NDC-classes-full": "#4daf4a",
        "NDC-substances-full": "#4daf4a",
        "tags-math-sx": "#984ea3",
        "tags-ask-ubuntu": "#984ea3",
        "threads-math-sx": "#ff7f00",
        "threads-ask-ubuntu": "#ff7f00",
        "coauth-MAG-Geology-full": "#FFC300",
        "coauth-MAG-History-full": "#FFC300",
        "coauth-DBLP-full": "#FFC300"
    }
    dataset_marker = {
        "email-Eu-full": "o",
        "email-Enron-full": "o",
        "contact-primary-school": "^",
        "contact-high-school": "^",
        "NDC-substances-full": "D",
        "NDC-classes-full": "D",
        "tags-math-sx": "P",
        "tags-ask-ubuntu": "P",
        "threads-ask-ubuntu": "s",
        "threads-math-sx": "s",
        "coauth-MAG-Geology-full": "<",
        "coauth-MAG-History-full": "<",
        "coauth-DBLP-full": "<"
    }
    for portion_str in portionlist:
        best_performing_alphas = []
        for idx, dataname in enumerate(dataset):
            best_performing_alpha = -1
            best_Dstat = 1000
            for algoname in algorithmlist:
                avg_degreeDstat = get_avg_degreeDstat(algoname, dataname, portion_str)
                if avg_degreeDstat < best_Dstat:
                    best_Dstat = avg_degreeDstat
                    alpha_str = algoname.split("_")[-1]
                    best_performing_alpha =  float(alpha_str)
            best_performing_alphas.append(best_performing_alpha)

        # Plot "Best-performing alpha - skewness" for all datasets
        plt.figure(figsize=(4.0, 3.1), dpi=120)
        xs, ys = [], []
        for didx, dname in enumerate(dataset):
            xs.append(math.log2(skewness[dname]))
            ys.append(math.log2(best_performing_alphas[didx] + 1))
        X, Y = np.array(xs).reshape(-1, 1), np.array(ys).reshape(-1, 1)
        cor = np.corrcoef(X.reshape(-1), Y.reshape(-1))[0,1]
        reg = LinearRegression().fit(X, Y)
        coef, intercept = reg.coef_, reg.intercept_
        X = np.arange(min(xs) - 0.5, max(xs) + 0.5).reshape(-1,1)
        plt.plot(X, reg.predict(X), alpha=0.5, color="black", linewidth=3, linestyle="dashed", zorder=1)
        
        for didx, dname in enumerate(dataset):
            plt.scatter(math.log2(skewness[dname]), math.log2(best_performing_alphas[didx] + 1), zorder=2, color=dataset_color[dname], marker=dataset_marker[dname], s=500, alpha=1.0)

        plt.ylabel(r"$log_{2}(\alpha^{*} + 1)$", fontsize=22)
        plt.xlabel(r"$log_{2}$" + "(skewness)", fontsize=22)
        
        ax = plt.gca()
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        plt.ylim(-1, ymax*1.1)
        ax.tick_params(labelcolor='#4B4B4B', labelsize=18)
        ax.text(0.4, 0.95, "corr.=%.2f" % (cor), transform=ax.transAxes, fontsize=18, verticalalignment='top')
        plt.tight_layout()
        portion = float(portion_str)
        plt.tight_layout()
        savedir = "figures/Observation/"
        if os.path.isdir(savedir) is False:
            os.makedirs(savedir)
        savefname = savedir + "skewness_%.1f.jpg" % (float(portion_str))
        plt.savefig(savefname, bbox_inches='tight')
        plt.close()

        
###################### Function ######################
def midas_observation3(dataname):

    ################### STYLE ######################
    dataset_color = {
        "email-Enron-full": "#e41a1c",
        "email-Eu-full": "#e41a1c",
        "contact-high-school": "#377eb8",
        "contact-primary-school": "#377eb8",
        "NDC-classes-full": "#4daf4a",
        "NDC-substances-full": "#4daf4a",
        "tags-math-sx": "#984ea3",
        "tags-ask-ubuntu": "#984ea3",
        "threads-math-sx": "#ff7f00",
        "threads-ask-ubuntu": "#ff7f00",
        "coauth-MAG-Geology-full": "#FFC300",
        "coauth-MAG-History-full": "#FFC300",
        "coauth-DBLP-full": "#FFC300"
    }
    dataset_marker = {
        "email-Eu-full": "o",
        "email-Enron-full": "o",
        "contact-primary-school": "^",
        "contact-high-school": "^",
        "NDC-substances-full": "D",
        "NDC-classes-full": "D",
        "tags-math-sx": "P",
        "tags-ask-ubuntu": "P",
        "threads-ask-ubuntu": "s",
        "threads-math-sx": "s",
        "coauth-MAG-Geology-full": "<",
        "coauth-MAG-History-full": "<",
        "coauth-DBLP-full": "<"
    }
    ################################################
    algorithmlist = ["es/global_deg_min_%.4f" % (2 ** i) for i in np.arange(-3,6.5,0.5)]
    algorithmlist += ["es/global_deg_min_0.0000"]

    best_performing_alphas = []
    for portion_str in portionlist:
        best_performing_alpha = -1
        best_Dstat = 1000
        for algoname in algorithmlist:
            avg_degreeDstat = get_avg_degreeDstat(algoname, dataname, portion_str)
            if avg_degreeDstat < best_Dstat:
                best_Dstat = avg_degreeDstat
                alpha_str = algoname.split("_")[-1]
                best_performing_alpha =  float(alpha_str)
        best_performing_alphas.append(best_performing_alpha)

    # Plot "best-performing alpha - portion" for each dataset
    xs, ys = [], []
    for pidx, portion_str in enumerate(portionlist):
        ys.append(math.log2(best_performing_alphas[pidx]+1))
        xs.append(int(100 * float(portion_str)))
    if dataname == "coauth-MAG-Geology-full":
        plt.figure(figsize=(4.3, 3.2), dpi=120)
    else:
        plt.figure(figsize=(4.2, 3.2), dpi=120)

    X, Y = np.array(xs).reshape(-1, 1), np.array(ys).reshape(-1, 1)
    cor = np.corrcoef(X.reshape(-1), Y.reshape(-1))[0,1]
    reg = LinearRegression().fit(X, Y)
    coef, intercept = reg.coef_, reg.intercept_
    xmin, xmax = min(xs), max(xs)
    X = np.arange(xmin * 0.9, xmax * 1.05).reshape(-1,1)
    plt.plot(X, reg.predict(X), alpha=0.5, color="black", linewidth=3, linestyle="dashed", zorder=1)

    plt.scatter(xs,  ys, marker=dataset_marker[dataname], alpha=1.0, s=550, color=dataset_color[dataname], zorder=2)
    plt.xlabel("Sampling Portion (%)", fontsize=22)
    plt.ylabel(r"$log_{2}(\alpha^{*} + 1)$", fontsize=22)

    ax = plt.gca()
    ax.tick_params(labelcolor='#4B4B4B', labelsize=18)
    ymin, ymax = ax.get_ylim()
    plt.ylim(ymin*0.9, ymax*1.1)
    ax.set_xticks([10, 20, 30, 40, 50])
    ax.text(0.4, 0.95, "corr.=%.2f" % (cor), transform=ax.transAxes, fontsize=18, verticalalignment='top')
    plt.tight_layout()
    savedir = "figures/Observation/alpha_portion/"
    if os.path.isdir(savedir) is False:
        os.makedirs(savedir)
    savefname = savedir + dataname + ".jpg"
    plt.savefig(savefname, bbox_inches='tight')
    plt.show()
    plt.close()

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--type', required=False, default='compete', type=str, help='baselines/compete/ablation')
    parser.add_argument('--select', required=True, type=int, help='0: Correlation Matrix, 1: Correlation With Degree, 2: MiDaS Observation1, 3: MiDaS Observation2, 4: MiDaS Observation3')
    parser.add_argument('--dataname', required=False, default='email-Eu-full', type=str, help="data name")
    parser.add_argument('--portion', required=False, default=0.3, type=float, help='[0.1,0.2,0.3,0.4,0.5]')
    args = parser.parse_args()
    
    if args.select == 0:
        correlation_matrix()
    elif args.select == 1:
        correlation_with_degree(args.portion)
    elif args.select == 2:
        midas_observation1(args.dataname, args.portion)
    elif args.select == 3:
        midas_observation2()
    elif args.select == 4:
        midas_observation3(args.dataname)
    