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
import matplotlib as mpl
plt.rcParams.update({'font.size': 13})

dataset = ["email-Enron-full", "email-Eu-full",
           "contact-high-school", "contact-primary-school",
           "NDC-classes-full", "NDC-substances-full",
           "tags-ask-ubuntu", "tags-math-sx",
           "threads-ask-ubuntu", "coauth-MAG-History-full", "coauth-MAG-Geology-full"]

algorithmlist=["ns/global_deg_0.0000", "ns/global_deg_1.0000", 
               "rw/rw_c_1", "ff/ff_c_0.51_0.20",
               "es/global_deg_min_0.0000", "tihs",
               "midas", "midas_grid",
               "mgs/add_degree", "mgs/add_avg", "mgs/exchange_degree", 
               "mgs/exchange_avg","mgs/remove_degree", "mgs/remove_avg"]

algorithm2labeling = {
    "ns/global_deg_0.0000": "RNS",
    "ns/global_deg_1.0000": "DNS",
    "rw/rw_c_1": "RW",
    "ff/ff_c_0.51_0.20": "FF",
    "es/global_deg_min_0.0000": "RHS",
    "tihs": "TIHS",
    "midas": "MiDaS",
    "midas_grid": "MiDaS-Grid",
    
    "mgs/add_degree": "MGS-DA",
    "mgs/add_avg": "MGS-AA",
    "mgs/exchange_degree": "MGS-DR",
    "mgs/exchange_avg": "MGS-AR",
    "mgs/remove_degree":"MGS-DD",
    "mgs/remove_avg": "MGS-AD",
}

def get_time(algoname, dataname, portion_str):
    fname = "../results_time/" + algoname + "/" + dataname + "_" + portion_str + "/time.txt"
    if os.path.isfile(fname) is False:
        print(fname)
        return -1
    else:
        with open(fname, "r") as f:
            _time = f.readline()
            agg_time = float(_time)

        return agg_time

if __name__ == "__main__":
    algo2timelist = defaultdict(float)
    for algoname in algorithmlist:
        for portion_str in ["0.10","0.20", "0.30", "0.40", "0.50"]: #["0.30"]:
            for data in dataset:
                time = get_time(algoname, data, portion_str)
                if time == -1 or algo2timelist[algoname] == -1:
                    print(algoname, data, portion_str)
                    continue
                else:
                    algo2timelist[algoname] += (time)
    print(algo2timelist)
    
    print(algo2timelist["midas_grid"] / algo2timelist["midas"])
    print(algo2timelist["ff/ff_c_0.51_0.20"] / algo2timelist["midas"])
    print(algo2timelist["mgs/remove_avg"] / algo2timelist["midas"])
    
    plt.figure(figsize=(8.0,4.4), dpi = 120)
    algonames = sorted(algorithmlist, key=lambda x: algo2timelist[x])
    for algoname in algonames:
        if "midas" == algoname:
            plt.bar(algorithm2labeling[algoname], algo2timelist[algoname], color="#be32be") #"#e41a1c")
        else:
            plt.bar(algorithm2labeling[algoname], algo2timelist[algoname], color="#b4b4b4") #"#377eb8")
    # plt.yscale('log', base=2)
    plt.yscale('log')
    plt.ylabel("Sampling Time\n(millisec.)", fontsize=20)
    ax = plt.gca()
    ax.tick_params (axis = 'x', labelrotation =72, labelsize=16)
    ax.tick_params('y', labelcolor='#4B4B4B', labelsize=16)

    my_colors = []
    for algoname in algonames:
        if "midas" == algoname:
            my_colors.append("#be32be")
        else:
            my_colors.append("black")
    for ticklabel, tickcolor in zip(plt.gca().get_xticklabels(), my_colors):
        ticklabel.set_color(tickcolor)

    minx, maxx = ax.get_xlim()
    maxy = max(algo2timelist.values())
    plt.hlines(y = maxy, xmin = algorithm2labeling["ff/ff_c_0.51_0.20"], xmax = maxx,
                   color = "#be32be", linestyle = 'dotted')
    plt.hlines(y = algo2timelist["midas"], xmin = algorithm2labeling["ff/ff_c_0.51_0.20"], xmax = maxx,
                   color = "#be32be", linestyle = 'dotted')

    plt.hlines(y = algo2timelist["midas"], xmin = algorithm2labeling["mgs/exchange_degree"], xmax = algorithm2labeling["mgs/add_degree"],
                   color = "#be32be", linestyle = 'dotted')
    plt.hlines(y = algo2timelist["midas_grid"], xmin = algorithm2labeling["mgs/exchange_degree"], xmax = algorithm2labeling["mgs/add_degree"],
                   color = "#be32be", linestyle = 'dotted')

    d = .015 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    # kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    # ax.plot((0-d,0+d), (1-(10*d),+1-(8*d)), **kwargs)
    # ax.plot((0-d,0+d),(1-(12*d),+1-(10*d)), **kwargs)
    # ax.set_yticks([2**18, 2**21, 2**24, 2**27])

    ax.annotate("",
                xy=(algorithm2labeling["midas_grid"], algo2timelist["midas"]*0.75), xycoords='data',
                xytext=(algorithm2labeling["midas_grid"], algo2timelist["midas_grid"]*1.25), textcoords='data',
                arrowprops=dict(arrowstyle="<->",
                                connectionstyle="arc3", color="#be32be", lw=2),
                )
    # plt.text(algorithm2labeling["rw/rw_c_1"], algo2timelist["mh/remove_avg"]*0.1, '     ' + r'$\ll$' + "\n  " + r'$\times \bf{839}$', color="#be32be", fontsize = 16)
    plt.text(algorithm2labeling["rw/rw_c_1"], algo2timelist["mgs/remove_avg"]*0.16,'   ' + r'$\times \bf{1018}$', color="#be32be", fontsize = 18)

    ax.annotate("",
                xy=(algorithm2labeling["mgs/remove_avg"], algo2timelist["mgs/remove_avg"]*1.2), xycoords='data',
                xytext=(algorithm2labeling["mgs/remove_avg"], algo2timelist["midas"]*0.75), textcoords='data',
                arrowprops=dict(arrowstyle="<->",
                                connectionstyle="arc3", color="#be32be", lw=2),
                )
    plt.text(algorithm2labeling["midas"], algo2timelist["midas_grid"] * 1.5, '     ' + r'$\times \bf{2.7}$', color="#be32be", fontsize = 18)

    ax.annotate("",
                xy=(algorithm2labeling["midas"], algo2timelist["midas"]*25), xycoords='data',
                xytext=(algorithm2labeling["midas"], algo2timelist["midas"]*0.8), textcoords='data',
                arrowprops=dict(arrowstyle="<-",
                                connectionstyle="arc3", color="#be32be", lw=2),
                )
    plt.text(algorithm2labeling["ns/global_deg_1.0000"], algo2timelist["midas"] * 32, '    ' + r'$\bf{MiDaS}$', color="#be32be", fontsize = 22)

    #     plt.legend(bbox_to_anchor=(1.05, 1))
    #     plt.title(portion_str)
    plt.tight_layout()
    savedir = "figures/Compete/"
    if os.path.isdir(savedir) is False:
        os.makedirs(savedir)
    savename = savedir + "sum_time.jpg"
    plt.savefig(savename, bbox_inches='tight')
    print(savename)
    