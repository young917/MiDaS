import pandas as pd
import numpy as np
from collections import defaultdict
from scipy import stats
import matplotlib.pyplot as plt
import os
import math
import matplotlib as mpl
import matplotlib.ticker as ticker
plt.rcParams.update({'font.size': 20})

baseline_algorithmlist=["ns/global_deg_0.0000", "ns/global_deg_1.0000", 
               "rw/rw_c_1", "ff/ff_c_0.51_0.20",
              "es/global_deg_min_0.0000", "tihs"]
compete_algorithmlist=["ns/global_deg_0.0000", "ns/global_deg_1.0000",

               "rw/rw_c_1", "ff/ff_c_0.51_0.20",
               "es/global_deg_min_0.0000", "tihs", 
               "mh/add_degree", "mh/add_avg", "mh/exchange_degree", "mh/exchange_avg",
               "mh/remove_degree", "mh/remove_avg",
               "midas"]

dataset = ["email-Enron-full", "email-Eu-full",
           "contact-high-school", "contact-primary-school",
           "NDC-classes-full", "NDC-substances-full",
           "tags-ask-ubuntu", "tags-math-sx", "threads-ask-ubuntu",
          "coauth-MAG-Geology-full","coauth-MAG-History-full"]

portionlist = ["0.10", "0.20", "0.30", "0.40", "0.50"]

############################ STYLE ################################
color_dict = {
    "es_skewness_search_ext": "#be32be", #"#942894", # 파란색

    "es/global_deg_min_0.0000":  "#377eb8", # 하늘색
    "ns/global_deg_0.0000": "#e41a1c",
    "ns/global_deg_1.0000": "#FF7A85", #"#F56E6E",

    "tihs": "#4daf4a", # 초록색
    "rw/rw_c_1": "#ff7f00", # 주황색
    "ff/ff_c_0.51_0.20": "#FFC81E", # 노랑색

    "mh/add_degree": "#CD853F", # 베이지색
    "mh/add_avg" : "#D7A35D",
    "mh/remove_degree" : "#D2691E", # 황토색
    "mh/remove_avg" : "#D7873C",
    "mh/exchange_degree" : "#8B4513", # 진한 갈색
    "mh/exchange_avg" : "#8B6331",

    # ANSWER
    'answer': "#000000"
}
line_style_dict = {
    "es/global_deg_min_0.0000": "solid",
    "es_skewness_search_ext": "solid",
    "ns/global_deg_0.0000": "solid",
    "ns/global_deg_1.0000": "solid", #"#F56E6E",
    
    "tihs": "solid",
    "rw/rw_c_1": "solid",
    "ff/ff_c_0.51_0.20": "solid",
    
    "mh/add_degree": "solid",
    "mh/add_avg" : (0, (1, 1)),
    "mh/remove_degree" : "solid",
    "mh/remove_avg" : (0, (1, 1)),
    "mh/exchange_degree" : "solid",
    "mh/exchange_avg" : (0, (1, 1)),
    
    # ANSWER    
    'answer': "(0,(5,1))" # 검정색
}
marker_dict = {
    "es/global_deg_min_0.0000": "o",
    "es_skewness_search_ext": "o",
    "ns/global_deg_0.0000": "o",
    "ns/global_deg_1.0000": "o", #"#F56E6E",

    "tihs": "o",
    "rw/rw_c_1": "o",
    "ff/ff_c_0.51_0.20": "o",

    "mh/add_degree": "o",
    "mh/add_avg" : "o",
    "mh/remove_degree" : "o",
    "mh/remove_avg" : "o",
    "mh/exchange_degree" : "o",
    "mh/exchange_avg" : "o",

    # ANSWER
    'answer': "," # 검정색
}

########################## Function ##############################
def robustness():
	algo2portion_degdstat = defaultdict(list)
	algo2portion_ranklist = defaultdict(list)
	algo2portion_normlist = defaultdict(list)
	
	# Read Results
	for portion_str in portionlist:
		# Rank
		d = pd.read_csv("./csvs/" + portion_str + "_compete_agg.csv")
		for i, row in d.iterrows():
			algoname = row["algorithm"]
			if row["algo opt"] != "-":
				algoname += "/" + row["algo opt"]
			if row["eval opt"] != "-":
				algoname += "/" + row["eval opt"]
			algo2portion_degdstat[algoname].append(float(row["degree avg"]))
			algo2portion_ranklist[algoname].append(float(row["aggregate rank"]))
		
		# Norm
		d = pd.read_csv("./csvs/" + portion_str + "_compete_norm_agg.csv")
		for i, row in d.iterrows():
			algoname = row["algorithm"]
			if row["algo opt"] != "-":
				algoname += "/" + row["algo opt"]
			if row["eval opt"] != "-":
				algoname += "/" + row["eval opt"]
			algo2portion_normlist[algoname].append(float(row["avg"]))

	portions = [10, 20, 30, 40, 50]
	# Plot D-Statistics in Degree Distribution
	plt.figure(figsize=(4.7,3.9), dpi=120)
	for algoname in algo2portion_degdstat.keys():
		if "midas" in algoname:
			plt.plot(portions, algo2portion_degdstat[algoname], label=algoname, alpha=1.0, color=color_dict[algoname], linewidth=10, linestyle=line_style_dict[algoname], marker=marker_dict[algoname],  markersize=10)
		else:
			plt.plot(portions, algo2portion_degdstat[algoname], label=algoname, alpha=1.0, color=color_dict[algoname], linewidth=3, linestyle=line_style_dict[algoname], marker=marker_dict[algoname],  markersize=5) # \ c= "#FF5675", \
	#plt.legend(fontsize=10, bbox_to_anchor=(-0.3, 1))
	ax = plt.gca()
	ax.tick_params(labelcolor='#4B4B4B', labelsize=18)
	ax.yaxis.set_label_coords(-0.22,0.5)
	ax.set_xticks(portions)

	plt.ylabel("D-Statistics in\nDegree Distribution", fontsize=20)
	plt.xlabel("Sampling Portion (%)", fontsize=20)
	plt.tight_layout()
	savedir = "figures/Compete/"
	if os.path.isdir(savedir) is False:
		os.makedirs(savedir)
	savename = savedir + "degree_dstat_portion.jpg"
	plt.savefig(savename, bbox_inches='tight')
	plt.close()

	# Plot Rank
	plt.figure(figsize=(4.0,3.9), dpi=120)
	for algo in algo2portion_ranklist.keys():
		if "midas" in algo:
			plt.plot(portions, algo2portion_ranklist[algo], marker=marker_dict[algo], markersize=10, label=algo, c=color_dict[algo],alpha=1.0, linewidth=10, linestyle=line_style_dict[algo])
		else:
			plt.plot(portions, algo2portion_ranklist[algo], marker=marker_dict[algo], markersize=5, label=algo, c=color_dict[algo], alpha=1.0, linewidth=3, linestyle=line_style_dict[algo])

	plt.ylabel("Avg. Rank", fontsize=20)
	plt.xlabel("Sampling Portion (%)", fontsize=20)
	ax = plt.gca()
	ax.tick_params(labelcolor='#4B4B4B', labelsize=18)
	ax.set_xticks(portions)
	# plt.legend(bbox_to_anchor=(-0.2, 1), fontsize=15)
	plt.tight_layout()
	savedir = "figures/Compete/"
	if os.path.isdir(savedir) is False:
		os.makedirs(savedir)
	savefname = savedir + "avgranking_portion.jpg"
	plt.savefig(savefname, bbox_inches='tight')
	plt.close()
	
	# Plot Norm
	plt.figure(figsize=(4.4,3.9), dpi=120)
	for algo in algo2portion_normlist.keys():
		if "midas" in algo:
			plt.plot(portions, algo2portion_normlist[algo], marker=marker_dict[algo],  markersize=10, label=algo, c=color_dict[algo],alpha=1.0, linewidth=10, linestyle=line_style_dict[algo])
		else:
			plt.plot(portions, algo2portion_normlist[algo], marker=marker_dict[algo],   markersize=5, label=algo, c=color_dict[algo], alpha=1.0, linewidth=3, linestyle=line_style_dict[algo])

	plt.ylabel("Avg. Z-Score", fontsize=20)
	plt.xlabel("Sampling Portion (%)", fontsize=20)
	ax = plt.gca()
	ax.tick_params(labelcolor='#4B4B4B', labelsize=18)
	# plt.legend(bbox_to_anchor=(-0.2, 1), fontsize=15)
	ax.set_xticks(portions)
	ax.yaxis.set_label_coords(-0.2,0.5)

	plt.tight_layout()
	savedir = "figures/Compete/"
	if os.path.isdir(savedir) is False:
		os.makedirs(savedir)
	savefname = savedir + "avgz_portion.jpg"
	plt.savefig(savefname, bbox_inches='tight')


########################## Function ##############################
labeldictx = {
    "degree": "Degree",
    "intersect" : "Int. Size",
    "pairdeg": "Pair degree",
    "size": "Size",
    "singular_value" : "Component",
    "size_wcc": "Rank"
}
labeldicty = {
    "degree": "Cumulative\nProbability",
    "intersect" :"Cumulative\nProbability",
    "pairdeg": "Cumulative\nProbability",
    "size": "Cumulative\nProbability",
    "singular_value" : "Relative\nVariance",
    "size_wcc": "# Nodes\n ",
    "density": "Density\n ",
    "overlapness": "Overlapness\n ",
    "effective_diameter" : "Diameter\n ",
    "global_cc": "GCC\n "
}
evallist = ["degree", "intersect", "pairdeg", "size",
            "singular_value", "size_wcc",
            "density", "overlapness", "effective_diameter", "global_cc"]

def get_dirpath(dir_path):
    if "answer" in dir_path:
        return dir_path
    
    _dir_path = dir_path
    assert os.path.isfile(dir_path + "/agg_result_arg.txt"), dir_path
    with open(dir_path + "/agg_result_arg.txt", "r") as f:
        for line in f.readlines():
            ename, arg = line[:-1].split(" : ")
            if ename == "degree":
                _dir_path = dir_path + "/" + arg + "/"

    return _dir_path

def get_dist(dir_path):
    dir_path = get_dirpath(dir_path)
    its = {}
    _its = pd.read_csv(dir_path + "intersect_dist.txt").sort_values(by=['intersect'])
    its['value'] = list(_its['value'].cumsum())
    its['intersect'] = list(_its['intersect'])
    pairdeg = {}
    _pairdeg = pd.read_csv(dir_path + "pairdeg_dist.txt").sort_values(by=['pairdeg'])
    pairdeg['value'] = list(_pairdeg['value'].cumsum())
    pairdeg['pairdeg'] = list(_pairdeg['pairdeg'])
    size = {}
    _size = pd.read_csv(dir_path + "size_dist.txt").sort_values(by=['size'])
    size['value'] = list(_size['value'].cumsum())
    size['size'] = list(_size['size'])
    deg = {}
    _deg = pd.read_csv(dir_path +  "degree_dist.txt").sort_values(by=['degree'])
    deg['value'] = list(_deg['value'].cumsum())
    deg['degree'] = list(_deg['degree'])

    singular_value_dist  = {}
    tmp_dict = {}
    with open(dir_path + "singular_value_dist2.txt", "r") as f:
        for line in f.readlines():
            line = line[:-1]
            sv, var = line.split(" ")
            tmp_dict[float(sv)] = float(var)
    keylist = sorted(list(tmp_dict.keys()))
    values = []
    for k in keylist:
        values.append(tmp_dict[k])
    singular_value_dist["singular_value"] = keylist
    singular_value_dist["value"] = values

    size_wcc = {}
    _size_wcc = pd.read_csv(dir_path + "sizewcc.txt").sort_values(by=['size_wcc'])
    size_wcc['value'] = list(_size_wcc['size_wcc'] / _size_wcc['num_nodes'])
    assert size_wcc["value"][-1] == 1
    size_wcc['size_wcc'] = list(range(1, len(_size_wcc)+1))

    density_pd = pd.read_csv(dir_path + "density.txt")
    density = float(density_pd["num edges"] / density_pd["num nodes"])

    global_cc = -1
    with open(dir_path + "global_cc.txt", "r") as f:
        global_cc = f.readline()
        global_cc = float(global_cc)

    overlapness = -1
    with open(dir_path + "overlapness.txt", "r") as f:
        overlapness = f.readline()
        overlapness = float(overlapness)

    effective_diam = -1
    with open(dir_path + "effdiameter.txt", "r") as f:
        effective_diam = f.readline()
        effective_diam = float(effective_diam)

    dic = {'intersect': its, 'pairdeg': pairdeg, 'size': size, 'degree': deg,
            "singular_value" : singular_value_dist,
           "size_wcc": size_wcc,"global_cc" : global_cc, "density" : density,  "overlapness" : overlapness, "effective_diameter" : effective_diam}
    return dic

def plot_figure(dataname, algo_dist, portion):
    for eval_name in evallist 
        if eval_name in ["global_cc", "density", "effective_diameter", "overlapness"]:
            plt.figure(figsize=(4.4,3.6), dpi=120)
			# plt.figure(figsize=(4.3, 2.4), dpi=120)
        elif eval_name in ["singular_value", "size_wcc"]:
            plt.figure(figsize=(3.5,3.6), dpi=120)
			# plt.figure(figsize=(4.3, 2.7), dpi=120)
        else:
            plt.figure(figsize=(4.3,3.6), dpi=120)
			# plt.figure(figsize=(4.3, 2.7), dpi=120)

        for algo_name in algo_dist.keys():
            dist = algo_dist[algo_name]
            assert dist[eval_name] is not None, dataname + " / " + algo_name + " / " + eval_name
            if algo_name == "answer":
                if isinstance(dist[eval_name], float):
                    continue
                elif len(dist[eval_name][eval_name]) == 1:
                    plt.scatter(dist[eval_name][eval_name][0], dist[eval_name]['value'][0], s=300, label=algo_name, alpha=1.0, c=color_dict[algo_name])
                else:
                    plt.plot(dist[eval_name][eval_name], dist[eval_name]['value'], linewidth=12, linestyle=line_style_dict[algo_name], label=algo_name, alpha=1.0, c=color_dict[algo_name]) #, marker=marker_dict[algo_name], markersize=4)
            else:
                if isinstance(dist[eval_name], float):
                    plt.bar(labelalgo[algo_name], dist[eval_name], label=algo_name, alpha=1.0, color=color_dict[algo_name], align='center', width=1)
                elif len(dist[eval_name][eval_name]) == 1:
                    plt.scatter(dist[eval_name][eval_name][0], dist[eval_name]['value'][0], s=300, label=algo_name, alpha=1.0, c=color_dict[algo_name])
                else:
                    if algo_name == "midas":
                        plt.plot(dist[eval_name][eval_name], dist[eval_name]['value'], linewidth=6, linestyle=line_style_dict[algo_name], label=algo_name, alpha=1.0, c=color_dict[algo_name]) #, marker=marker_dict[algo_name], markersize=2.5)
                    else:
                        plt.plot(dist[eval_name][eval_name], dist[eval_name]['value'], linewidth=4, linestyle=line_style_dict[algo_name], label=algo_name, alpha=1.0, c=color_dict[algo_name]) #, marker=marker_dict[algo_name], markersize=2.5)
        ax = plt.gca()
        xmin, xmax = ax.get_xlim()
		if eval_name in ["degree", "intersect", "pairdeg", "size"]:
            plt.xscale('log')
			locmaj = mpl.ticker.LogLocator(numticks=5)
            ax.xaxis.set_major_locator(locmaj)
        if eval_name in ["global_cc", "density", "effective_diameter", "overlapness"]:
            plt.hlines(y=algo_dist["answer"][eval_name], xmin=xmin-0.02, xmax=xmax+0.02, linewidth=7, linestyle="dashdot", color='black', alpha=1.0)
            ymin, ymax = ax.get_ylim()
            plt.ylim(ymin, ymax*1.1)
            ax.set_xticks([])
			plt.xlabel( "Algorithm", fontsize=21)
            # ax.tick_params(axis='y',labelcolor='#4B4B4B', labelsize=19)
            # ax.tick_params(axis='x',labelcolor='white', labelsize=19)
		else:
            plt.xlabel(labeldictx[eval_name], fontsize=20)
        ax.tick_params(axis='y', which='major', labelcolor='#4B4B4B', labelsize=16)
        ax.tick_params(axis='x', which='minor', labelcolor='#4B4B4B',labelsize=14)
        ax.tick_params(axis='x', which='major', labelcolor='#4B4B4B', labelsize=16)

        plt.ylabel(labeldicty[eval_name], fontsize=21)
		# plt.legend(title= dataname, bbox_to_anchor=(-0., 1), loc='upper left')
        # if eval_name == "size_wcc":
            # ax.yaxis.set_label_coords(-0.38,0.5)
        # elif eval_name in ["global_cc", "density", "effective_diameter", "overlapness"]:
            # ax.yaxis.set_label_coords(-0.38,0.5)
        # else:
            # ax.yaxis.set_label_coords(-0.27,0.5)

		plt.tight_layout()
		if type == "baseline":
			savedir = "figures/Baseline/" + dataname + "/%.1f"% (portion) + "/"
			savefname = "figures/Baseline/" + dataname + "/%.1f"% (portion) + "/dist_" + eval_name + ".jpg"
		else:
        	savedir = "figures/Compete/" + "%.1f"% (portion) + "/"
        	savefname = savedir + "dist_" + eval_name + "_" + dataname + ".jpg"
		if os.path.isdir(savedir) is False:
            os.makedirs(savedir)
        plt.savefig(savefname, bbox_inches='tight')
        plt.show()
        plt.close()

def plot_dist(type="compete"):
	if type == "baseline":
		algorithmlist = baseline_algorithmlist
	else:
		algorithmlist = compete_algorithmlist

	for portion_str in portionlist:
		for dataname in dataset:
			algo2dist = {}
			algo2dist["answer"] = get_dist("../results/answer_dist/" + dataname + "/")
			for algo_name in algorithmlist:
				dir_path = "../results/" + algo_name
				algo_dir = dir_path + "/" + dataname + "_" + portion_str + "/"
				ret = get_dist(algo_dir)
				if ret == -1:
					continue
				algo2dist[algo_name] = ret
			print(dataname)
			plot_figure(dataname, algo2dist, float(portion_str), type)


