import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
from itertools import chain
import os
import argparse
import math
from collections import defaultdict
import shutil
from copy import deepcopy
from scipy.stats import norm
import math
from sklearn.linear_model import LinearRegression

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
    "coauth-MAG-History-full" : 52.414521,
    "coauth-MAG-Geology-full" : 12.639768,
}

dataset = ["contact-high-school", "contact-primary-school", "email-Enron-full", "email-Eu-full", "tags-math-sx", "tags-ask-ubuntu", "NDC-classes-full", "NDC-substances-full", "coauth-MAG-Geology-full", "coauth-MAG-History-full", "threads-ask-ubuntu"]

# portionlist = ["0.10","0.20","0.30", "0.4s0", "0.50"]
portionlist = [0.1, 0.2, 0.3, 0.4, 0.5]

def _read(fname):
    deg = -1
    if os.path.isfile(fname) is False:
        return -1
    with open(fname , "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.split("\n")[0]
            tmp = line.split(" : ")
            if len(tmp) == 2:
                eval_name, result = tmp
                if eval_name == "degree":
                    deg = float(result)
                    break
    return deg

def query_degDstat_Time(algoname, data, portion):
    portion_str = "%.2f" % (portion)
    dir_path = "../results/" + algoname + "/" + data + "_" + portion_str
    avg_degDstat = 0.0
    avg_time = 0.0
    count = 0
    for i in range(1,4):
        fname = dir_path + "/" + str(i) + "/result.txt"
        degDstat = _read(fname)
        if degDstat == -1:
#             print(fname)
            continue
        avg_degDstat += degDstat
        count += 1
    if count == 0:
        print(algoname + "/" + data + "_" + portion_str)
        return -1, -1

    avg_degDstat /= count

    fname = "../results_time/" + algoname + "/" + data + "_" + portion_str + "/adjust_time.txt"
    with open(fname, "r") as f:
        avg_time = float(f.readline())

    return avg_degDstat, avg_time

def get_alpha(algoname):
    alpha_str = algoname.split("_")[-1]
    return float(alpha_str)

def train_linear_regressor(except_dataset, grid_result_alpha):
    X, Y = [], []
    except_domain = except_dataset.split("-")[0]
    for portion in portionlist:    
        for idx, dname in enumerate(dataset):
            #if dname == except_dataset:
            if except_domain in dname:
                continue
            X.append([portion, math.log2(skewness[dname]) ])
            Y.append(math.log2(grid_result_alpha[portion][dname] + 1))

    X, Y = np.array(X), np.array(Y)
    reg = LinearRegression().fit(X, Y)
    score = reg.score(X, Y)
    coef, intercept = reg.coef_, reg.intercept_
    
    return reg, coef, intercept

def trained_model(reg, portion, sk):
    x = np.array( [[portion, math.log2(sk)]] )
    predicted_logscale= reg.predict(x)
    predicted = 2 ** predicted_logscale - 1
    
    return predicted_logscale, predicted

def get_fname(alpha):
    return "es/global_deg_min_%.4f" % (alpha)
    
def run_midas_grid(portion, data, search_space):
    algorithmlist = ["es/global_deg_min_%.4f" % a for a in search_space]

    best_performing_alpha = -1
    best_dstat = 1000
    time = 0
    trynum = 0
    for algoname in algorithmlist:
        degDstat, time = query_degDstat_Time(algoname, data, portion)
        time += time
        trynum += 1
        if degDstat < best_dstat:
            best_dstat = degDstat
            best_performing_alpha = get_alpha(algoname)
            
    # Save Time
    portion_str = "%.2f" % (portion)
    outputdir = "../results_time/midas_grid/" + data + "_" + portion_str + "/"
    if os.path.isdir(outputdir) is False:
        os.makedirs(outputdir)
    with open(outputdir + "time.txt", "w") as f:
        f.write(str(time) + "\n")

    return best_performing_alpha, best_dstat, time, trynum

def run_midas(data, portion, trained_linear_model, search_space):
    best_performing_alpha, best_dstat, time, trynum = -1, 1000, 0, 0
    _, _alpha = trained_model(trained_linear_model, portion, skewness[data])

    start_idx = -1
    while start_idx < (len(search_space) - 1):
        if search_space[start_idx + 1] < _alpha:
            start_idx += 1
        else:
            break

    if start_idx == -1:
        search_direction = 1
        start_idx = 0
    elif start_idx == len(search_space) - 1:
        search_direction = -1
    else:
        left_alpha = search_space[start_idx]
        degDstat_left_alpha, spent = query_degDstat_Time(get_fname(left_alpha), data, portion)
        time += spent
        trynum += 1

        right_alpha = search_space[start_idx + 1]
        degDstat_right_alpha, spent = query_degDstat_Time(get_fname(right_alpha), data, portion)
        time += spent
        trynum += 1

        if degDstat_left_alpha <= degDstat_right_alpha:
            best_performing_alpha = left_alpha
            best_dstat = degDstat_left_alpha
            search_direction = -1
            start_idx = start_idx - 1
        else:
            best_performing_alpha = right_alpha
            best_dstat = degDstat_right_alpha
            search_direction = 1
            start_idx = start_idx + 2

    idx = start_idx
    while (idx >= 0) and (idx < len(search_space)):
        next_alpha = search_space[idx]
        degDstat, spent = query_degDstat_Time(get_fname(next_alpha), data, portion)
        time += spent
        trynum += 1
        if best_dstat <= degDstat:
            break
        best_performing_alpha = next_alpha
        best_dstat = degDstat
        idx += search_direction

    return best_performing_alpha, best_dstat, time, trynum

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', required=True, type=str)
    parser.add_argument('--portion', required=True, type=float)
    args = parser.parse_args()

    search_space = [float("%.4f" % (2 ** i)) for i in np.arange(-3,6.5,0.5)] + [0.0]
    
    # MiDaS_Grid
    result_alpha = defaultdict(dict)
    result_dstat = defaultdict(dict)
    result_time = defaultdict(dict)
    result_try = defaultdict(dict)
    for portion in portionlist:
        for idx, data in enumerate(dataset):
            best_performing_alpha, best_dstat, time, trynum = run_midas_grid(portion, data, search_space)
            result_alpha[portion][data] = best_performing_alpha
            result_dstat[portion][data] = best_dstat
            result_time[portion][data] = time
            result_try[portion][data] = trynum
    
    # print(result_alpha)
    midas_grid = {
        "alpha": result_alpha,
        "degDstat": result_dstat,
        "time": result_time,
        "try": result_try
    }
    # Save Result
    with open("../results/midas_grid/search_result.txt", "w") as f:
        for portion in portionlist:
            portion_str = "%.2f" % (portion)
            for data in dataset:
                best_performing_alpha = result_alpha[portion][data]
                alpha_name = "../results/es/global_deg_min_%.4f" % (best_performing_alpha)
                f.write(alpha_name + "/" + data + "_" + portion_str + "\n")

    # MiDaS
    trained_linear_model, coef, intercept = train_linear_regressor(args.data, midas_grid["alpha"])
    best_performing_alpha, best_dstat, time, trynum = run_midas(args.data, args.portion, trained_linear_model, search_space)

    outputdir = "../results_time/midas/"
    if os.path.isdir(outputdir) is False:
        os.makedirs(outputdir)
    with open(outputdir + "time.txt", "w") as f: # save time
        f.write(str(time) + "\n")

    outputdir = "../results/midas/"
    if os.path.isdir(outputdir) is False:
        os.makedirs(outputdir)
    outputname = outputdir + "search_result.txt"
    with open(outputname, "a+") as f: # save result
        f.write("../results/es/global_deg_min_%.4f" % (best_performing_alpha) + "/" + args.data + "_%.2f" % (args.portion) + "\n")

