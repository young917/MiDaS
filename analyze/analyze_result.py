import pandas as pd
import numpy as np
from collections import defaultdict
from scipy import stats
import matplotlib.pyplot as plt
import os
import math
import argparse

#################### Default  #######################################################

baseline_algorithmlist=["es/global_deg_min_0.0000", "tihs", "ns/global_deg_1.0000", "rw/rw_c_1", "ff/ff_c_0.51_0.20", "ns/global_deg_0.0000"]
compete_algorithmlist = ["es/global_deg_min_0.0000",  "tihs", "ff/ff_c_0.51_0.20", "ns/global_deg_0.0000", "ns/global_deg_1.0000", "mgs/add_degree", "mgs/add_avg", "mgs/exchange_degree", "mgs/exchange_avg", "mgs/remove_degree", "mgs/remove_avg", "rw/rw_c_1", "midas"]
ablation_algorithmlist=["midas_grid_ablation", "avgdegree", "maxdegree", "midas_ns"]

portionlist = ["0.10", "0.20", "0.30", "0.40", "0.50"]

dataset = ["email-Enron-full", "email-Eu-full", "contact-high-school", "contact-primary-school", "NDC-classes-full", "NDC-substances-full", "tags-ask-ubuntu", "tags-math-sx", "threads-ask-ubuntu", "coauth-MAG-Geology-full","coauth-MAG-History-full"]

evallist = ["degree", "intersect", "pairdeg", "size", "Time", "singular_value", "size_wcc", "global_cc_norm", "density_norm", "overlapness_norm", "effective_diameter_norm"]

def get_evaluation(dir_path):
    evaluation = {}
    inputfname = dir_path + "/agg_entire_evaluation.txt"
    if not os.path.isfile(inputfname):
        print("No " + inputfname)
        return -1
    
    with open(inputfname , "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.split("\n")[0]
            tmp = line.split(" : ")
            if len(tmp) == 2:
                eval_name, result = tmp
                if eval_name in evallist:
                    evaluation[eval_name] = float(result)
    for ename in evallist:
        assert ename in evaluation, ename + "  " + dir_path
            
    return evaluation

def aggregate_algorithms(portion_str, dataname, type):
    if type == "baseline":
        algorithmlist = baseline_algorithmlist
    elif type == "compete":
        algorithmlist = compete_algorithmlist
    elif type == "ablation":
        algorithmlist = ablation_algorithmlist
        
    algo2evaldict = {}
    # Make Output
    outputname = "./csvs/" + dataname + "/" + portion_str + "_" + type + ".csv"
    if os.path.isdir("./csvs/" + dataname + "/") is False:
        os.makedirs("./csvs/" + dataname + "/")
    with open(outputname, "w") as f:
        f.write("algorithm,algo opt,eval opt,")
        line = ",".join(evallist)
        f.write(line + "\n")
        
    # Save in algo2evaldict
    for algoname in algorithmlist:
        ret_eval = get_evaluation("../results/" + algoname + "/" + dataname + "_" + portion_str)
        algo2evaldict[algoname] = ret_eval
    
    for algoname in algo2evaldict.keys():
        # print(algoname)
        tmp = algoname.split("/")
        if len(tmp) == 3:
            algorithm, algo_opt, eval_opt = tmp
        elif len(tmp) == 2:
            algorithm, algo_opt = tmp
            eval_opt = "-"
        elif len(tmp) == 1:
            algorithm = tmp[0]
            algo_opt = "-"
            eval_opt = "-"
        evaldict = algo2evaldict[algoname]
        if evaldict == -1:
            print("No " + dataname + " " + algoname)
        else:
            with open(outputname, "a") as f:
                strlist = []
                for ename in evallist:
                    if math.fabs(evaldict[ename]) < 1:
                        strlist.append("%.18f" % (evaldict[ename]))
                    else:
                        strlist.append(str(evaldict[ename]))
                line = ",".join([algorithm, algo_opt, eval_opt] + strlist)
                f.write(line + "\n")
    
    # Make Normalize Result
    d = pd.read_csv("./csvs/" + dataname + "/" + portion_str + "_" + type + ".csv")
    for col in ["global_cc_norm", "overlapness_norm", "density_norm", "effective_diameter_norm"]:
        d[col] = d[col].abs()
    for col in list(d.columns)[3:]:
        if d[col].std() != 0:
            d[col] = (d[col] - d[col].mean()) / d[col].std()
    # CHANGED!!
    selected_columns = [e for e in evallist if e != 'Time']
    norms = d[selected_columns]
    d['avg'] = norms.mean(axis=1)
    d = d.sort_values(by=["avg"], ascending=True)
    d.to_csv("./csvs/" + dataname + "/" + portion_str + "_" + type + "_norm.csv", index=False)

    # Make Ranking Result
    d = pd.read_csv("./csvs/" + dataname + "/" + portion_str + "_" + type + ".csv")
    for ename in (evallist):
        d['rank ' + ename] = d[ename].abs().rank(method='min')
    selected_columns = ['rank ' + ename for ename in evallist if ename != "Time"]
    ranks = d[selected_columns]
    d['eval avg rank'] = ranks.mean(axis=1)
    d = d.sort_values(by=["eval avg rank"], ascending=True)
    d.to_csv("./csvs/" + dataname + "/" + portion_str + "_" + type + ".csv", index=False)


def aggregate_all_datasets(portion_str, type):
    ######  Normalize ######
    algoname2normevallist = {}
    algoname2count = defaultdict(int)    
    norm_columns = evallist
    for didx, dataname in enumerate(dataset):
        d = pd.read_csv("./csvs/" + dataname + "/" + portion_str + "_" + type + "_norm.csv")
        for i, row in d.iterrows():
            algoname = row["algorithm"]
            if row["algo opt"] != "-":
                algoname += "/" + row["algo opt"]
            if row["eval opt"] != "-":
                algoname += "/" + row["eval opt"]
            if algoname not in algoname2normevallist:
                algoname2normevallist[algoname] = defaultdict(list)
            for idx, col in enumerate(norm_columns):
                algoname2normevallist[algoname][col].append(float(row[col]))
            algoname2count[algoname] += 1

    for algoname in algoname2normevallist.keys():
        #print(algoname)
        #print(algoname2normevallist[algoname]["Time"], sum(algoname2normevallist[algoname]["Time"]) / len(algoname2normevallist[algoname]["Time"]))
        for col in norm_columns:
            algoname2normevallist[algoname][col] = np.mean(algoname2normevallist[algoname][col])
        #print(algoname2normevallist[algoname]["Time"])
    norm_outputname = "./csvs/" + portion_str + "_" + type + "_norm_agg.csv"
    with open(norm_outputname, "w") as f:
        f.write("algorithm,algo opt,eval opt,")
        line = ",".join(norm_columns)
        f.write(line + "\n")

    for algoname in algoname2normevallist.keys():
        # print(algoname)
        tmp = algoname.split("/")
        if len(tmp) == 3:
            algorithm, algo_opt, eval_opt = tmp
        elif len(tmp) == 2:
            algorithm, algo_opt = tmp
            eval_opt = "-"
        elif len(tmp) == 1:
            algorithm = tmp[0]
            eval_opt = "-"
            algo_opt = "-"
        #count = algoname2count[algoname]
        assert algoname2count[algoname] == len(dataset), algoname
        with open(norm_outputname, "a") as f:
            line = ",".join([algorithm, algo_opt, eval_opt] + [str(algoname2normevallist[algoname][e]) for e in norm_columns])
            f.write(line + "\n")
    
    d = pd.read_csv(norm_outputname)
    norm_target_cols = [e for e in evallist if e != 'Time']
    tmp = d[norm_target_cols]
    d["avg"] = tmp.mean(axis=1)
    d = d.sort_values(by=["avg"], ascending=True)
    d.to_csv(norm_outputname, index=False)

    # Ranking
    ranking_columns = evallist + ['rank ' + e for e in evallist]
    algoname2evallist = {}
    algoname2count = defaultdict(int)
    for didx, dataname in enumerate(dataset):
        # print(dataname)
        d = pd.read_csv("./csvs/" + dataname + "/" + portion_str + "_" + type + ".csv")
        for i, row in d.iterrows():
            algoname = row["algorithm"]
            if row["algo opt"] != "-":
                algoname += "/" + row["algo opt"]
            if row["eval opt"] != "-":
                algoname += "/" + row["eval opt"]
            if algoname not in algoname2evallist:
                algoname2evallist[algoname] = defaultdict(list)   
            # make a list for each feature in order of dataset
            for idx, col in enumerate(ranking_columns):
                algoname2evallist[algoname][col].append(math.fabs(float(row[col])))

            algoname2count[algoname] += 1

    for algoname in algoname2evallist.keys():
        for col in ranking_columns:
            algoname2evallist[algoname][col + " sd"] = np.std(algoname2evallist[algoname][col])
            algoname2evallist[algoname][col + " avg"] = np.mean(algoname2evallist[algoname][col])
    ranking_output_columns = []
    for ev in ranking_columns:
        if ev == "Time" or ev == "rank Time":
            ranking_output_columns.append(ev + " avg")
        else:
            ranking_output_columns.append(ev + " avg")
            ranking_output_columns.append(ev + " sd")
    ranking_outputname = "./csvs/" + portion_str + "_" + type + "_agg.csv"
    with open(ranking_outputname, "w") as f:
        line = ",".join(ranking_output_columns)
        f.write("algorithm,algo opt,eval opt," + line + "\n")
    
    for algoname in algoname2evallist.keys():
        tmp = algoname.split("/")
        if len(tmp) == 3:
            algorithm, algo_opt, eval_opt = tmp
        elif len(tmp) == 2:
            algorithm, algo_opt = tmp
            eval_opt = "-"
        elif len(tmp) == 1:
            algorithm = tmp[0]
            algo_opt = "-"
            eval_opt = "-"
            
        #count == algoname2count[algoname]
        assert algoname2count[algoname] == len(dataset), algoname + " " + str(algoname2count[algoname])
        
        with open(ranking_outputname, "a") as f:
            line = ",".join([algorithm, algo_opt, eval_opt] + [str(algoname2evallist[algoname][e]) for e in ranking_output_columns])
            f.write(line + "\n")
    
    d = pd.read_csv(ranking_outputname)
    ranking_target_cols = ['rank ' + e + ' avg' for e in evallist if e != 'Time']
    # ranking_target_cols = ['rank degree avg', 'rank size avg', 'rank intersect avg', 'rank pairdeg avg', 'rank size_wcc avg', 'rank global_cc_norm avg', 'rank singular_value avg' , 'rank density_norm avg', 'rank overlapness_norm avg', 'rank effective_diameter_norm avg']
    
    ranks = d[ranking_target_cols]
    d['aggregate rank'] = ranks.mean(axis=1)
    d = d.sort_values(by=["aggregate rank"], ascending=True)
    d.to_csv(ranking_outputname, index=False)

def aggregate_all_portions(type):
    algo2rank_prop = {}
    algo2norm_prop = {}
    algo2dstat_prop = {}

    for portion_str in portionlist:
        # Rank
        d = pd.read_csv("./csvs/" + portion_str + "_" + type + "_agg.csv")
        for i, row in d.iterrows():
            algoname = row["algorithm"]
            if row["algo opt"] != "-":
                algoname += "/" + row["algo opt"]
            if row["eval opt"] != "-":
                algoname += "/" + row["eval opt"]

            if algoname not in algo2rank_prop:
                algo2rank_prop[algoname] = defaultdict(list)
                algo2dstat_prop[algoname] = defaultdict(list)
            for evalname in evallist:
                algo2rank_prop[algoname][evalname].append(float(row["rank " + evalname + " avg"]))
                algo2dstat_prop[algoname][evalname].append(float(row[evalname + " avg"]))
            algo2rank_prop[algoname]["avg"].append(float(row["aggregate rank"]))
        # Norm
        d = pd.read_csv("./csvs/" + portion_str + "_" + type + "_norm_agg.csv")
        for i, row in d.iterrows():
            algoname = row["algorithm"]
            if row["algo opt"] != "-":
                algoname += "/" + row["algo opt"]
            if row["eval opt"] != "-":
                algoname += "/" + row["eval opt"]
                
            if algoname not in algo2norm_prop:
                algo2norm_prop[algoname] = defaultdict(list)
            for evalname in evallist:
                algo2norm_prop[algoname][evalname].append(float(row[evalname]))
            algo2norm_prop[algoname]["avg"].append(float(row["avg"]))
            
    for algo in algo2rank_prop:
        for evalname in evallist:
            algo2rank_prop[algo][evalname] = np.mean(algo2rank_prop[algo][evalname])
            algo2norm_prop[algo][evalname] = np.mean(algo2norm_prop[algo][evalname])
            algo2dstat_prop[algo][evalname] = np.mean(algo2dstat_prop[algo][evalname])
        algo2rank_prop[algo]['avg'] = np.mean(algo2rank_prop[algo]['avg'])
        algo2norm_prop[algo]['avg'] = np.mean(algo2norm_prop[algo]['avg'])

    # D-Stat Output
    dstat_outputname = "./csvs/allportion_" + type + "_agg_dstat.csv"
    columns = [e for e in evallist if e != 'Time']
    with open(dstat_outputname, "w") as f:
        line = ",".join(columns)
        f.write("algorithm,algo opt,eval opt,")
        f.write(line + "\n")
        
    for algoname in algo2dstat_prop.keys():
        tmp = algoname.split("/")
        if len(tmp) == 3:
            algorithm, algo_opt, eval_opt = tmp
        elif len(tmp) == 2:
            algorithm, algo_opt = tmp
            eval_opt = "-"
        elif len(tmp) == 1:
            algorithm = tmp[0]
            algo_opt = "-"
            eval_opt = "-"  
        else:
            print(tmp)
        evaldict = algo2dstat_prop[algoname]
        if evaldict == -1:
            print("Np " + dataname + " " + algoname)
        else:
            with open(dstat_outputname, "a") as f:
                line = ",".join([algorithm, algo_opt, eval_opt] + [ str(evaldict[ename]) for ename in evallist if ename != "Time"])
                f.write(line + "\n")
    
    # Rank Output
    ranking_outputname = "./csvs/allportion_" + type + "_agg_rank.csv"
    columns = [e for e in evallist if e != 'Time'] + ["avg"]
    with open(ranking_outputname, "w") as f:
        line = ",".join(columns)
        f.write("algorithm,algo opt,eval opt,")
        f.write(line + "\n")
        
    for algoname in algo2rank_prop.keys():
        tmp = algoname.split("/")
        if len(tmp) == 3:
            algorithm, algo_opt, eval_opt = tmp
        elif len(tmp) == 2:
            algorithm, algo_opt = tmp
            eval_opt = "-"
        elif len(tmp) == 1:
            algorithm = tmp[0]
            algo_opt = "-"
            eval_opt = "-"
        else:
            print(tmp)
        evaldict = algo2rank_prop[algoname]
        if evaldict == -1:
            print("Np " + dataname + " " + algoname)
        else:
            with open(ranking_outputname, "a") as f:
                line = ",".join([algorithm, algo_opt, eval_opt] + [ str(evaldict[ename]) for ename in columns if ename != "Time"])
                f.write(line + "\n")

    # Norm Output
    norm_outputname = "./csvs/allportion_" + type +"_agg_norm.csv"
    columns = [e for e in evallist if e != 'Time'] + ["avg"]
    with open(norm_outputname, "w") as f:
        line = ",".join(columns)
        f.write("algorithm,algo opt,eval opt,")
        f.write(line + "\n")
    for algoname in algo2norm_prop.keys():
        tmp = algoname.split("/")
        if len(tmp) == 3:
            algorithm, algo_opt, eval_opt = tmp
        elif len(tmp) == 2:
            algorithm, algo_opt = tmp
            eval_opt = "-"
        elif len(tmp) == 1:
            algorithm = tmp[0]
            algo_opt = "-"
            eval_opt = "-"
        else:
            print(tmp)
        evaldict = algo2rank_prop[algoname]
        # print(algoname)
        evaldict = algo2norm_prop[algoname]
        if evaldict == -1:
            print("Np " + dataname + " " + algoname)
        else:
            with open(norm_outputname, "a") as f:
                line = ",".join([algorithm, algo_opt, eval_opt]+ [ str(evaldict[ename]) for ename in evallist + ["avg"] if ename != "Time"])
                f.write(line + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--type', required=True, default='baseline', type=str, help='baselines/compete/ablation')
    args = parser.parse_args()

    for portion_str in portionlist:
        for dataname in dataset:
            aggregate_algorithms(portion_str, dataname, args.type)
        aggregate_all_datasets(portion_str, args.type)  
    aggregate_all_portions(args.type)

