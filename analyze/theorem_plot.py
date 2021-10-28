from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import argparse

dataset = ["email-Enron-full", "email-Eu-full", "contact-high-school", "contact-primary-school",
          "NDC-classes-full", "NDC-substances-full", "tags-ask-ubuntu", "tags-math-sx",
          "threads-ask-ubuntu", "threads-math-sx",
          "coauth-DBLP-full", "coauth-MAG-Geology-full", "coauth-MAG-History-full"]

color_red = "#e41a1c"
color_green = "#4daf4a" 
color_blue = "#377eb8"

def analyze(dname):
    degree_path = "../results/answer_dist/" + dname + "/degree_dist.txt"
    d = pd.read_csv(degree_path)
    degree_list = list(d["degree"]) 
    value_list = list(d["value"])
    degree_dict = {}
    avg_deg = 0
    for idx, deg in enumerate(degree_list):
        degree_dict[deg] = value_list[idx]
    for deg in sorted(degree_list):
        avg_deg += deg * degree_dict[deg]
        
    for opt_idx, opt in enumerate(["min", "max" ,"avg"]):
        print(dname, opt)
        plt.figure(figsize = (4.2,3), dpi=120)
        filepath = "../results/answer_dist/" + dname + "/phi_" + opt + "_analysis.txt"
        d = pd.read_csv(filepath)
        min_list, max_list, degree_list = list(d["min"]), list(d["max"]), list(d["degree"])
        final_mean = list(d["avg"])[-1]
        k_thres = min(degree_list) - 1
        for i, _max in enumerate(max_list):
            if _max <= final_mean:
                k_thres = max(k_thres, degree_list[i])
        plt.scatter(d["degree"], d["max"], s=100, c=color_blue)
        plt.hlines(final_mean, min(degree_list[0], k_thres), max(degree_list), linewidth=5, color=color_red)
        print(final_mean)
        
        ax = plt.gca()
        ax.tick_params(labelcolor='#4B4B4B', labelsize=16)
        if k_thres == min(degree_list) - 1:
            plt.vlines(min(degree_list) - 0.1, min_list[-1], max_list[-1], linewidth=5, linestyle="dashed", color=color_green)
            print("not exist")
        else:
            plt.vlines(k_thres, min_list[-1], max_list[-1], linewidth=5, color=color_green)
            print(k_thres)
            
        plt.vlines(avg_deg, min_list[-1], max_list[-1], linewidth=5, alpha=0.8, linestyle="dashed", color="gray")
        plt.xscale('log')
    
        if "contact" in dname:
            print(dname)
            locmaj = mpl.ticker.LogLocator(numticks=10)
            ax.xaxis.set_major_locator(locmaj)
        
        plt.xlabel("k", fontsize=20)
        plt.ylabel(r"$ln\phi(e)$", fontsize=20)
        plt.tight_layout()
        savedir = "figures/Theorem/" + opt + "/"
        if os.path.isdir(savedir) is False:
            os.makedirs(savedir)
        savefname = savedir + dname + ".jpg"
        plt.savefig(savefname, bbox_inches='tight')
        plt.show()
        plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', required=True, type=str, help="data name", default="email-Eu-full")
    args = parser.parse_args()
    
    analyze(args.data)
        