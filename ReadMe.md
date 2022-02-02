# MiDaS: Representative Sampling from Real-world Hypergraphs

We provide source code for (1) sampling hypergraphs, (2) evaluating the quality of sub-hypergraph and (3) reproducing all figures in the paper: MiDaS: Representative Sampling from Real-world Hypergraphs, Minyoung Choe, Jaemin Yoo, Geon lee, Woonsung Baek, U Kang and Kijung Shin, WWW 2022.

(1) **Sampling Hypergraphs**

* *Simple and Intuitive Approaches* : six intuitive approaches(RNS, RDN, RW, FF, RHS, TIHS) experimenting on 11 datasets with a focus on their limitations in preserving 10 properties
* *MiDaS-Basic* : a preliminary version of MiDaS which can improve the limitations of intuitive approaches by proper alpha setting
* *MiDaS*: a full-fledged version which can automatically tune alpha

(2) **Evaluating the Quality of Sub-Hypergraph**

Measure how precisely the sub-hypergraph preserves the structural properties of the entire hypergraph with respect to,

* *node-level statistics*: the distributions of node degrees and node-pair degrees
* *hyperedge-level statistics* : the distributions of hyperedge sizes and intersection sizes
* *graph-level statistics* : average clustering coefficient, density, overlapness, and effective diameter

(3) Reproducing *all* figures in the main paper

- - -

## Datasets

In the papaer, we used datasets after removing duplicated hyperedges. We preprocessed eleven datasets collected by [Austin R. Benson](https://www.cs.cornell.edu/~arb/data/). The datasets used in the paper are available in the "dataset" folder.

- - -

## How to Run

### Example

You can run all thirteen sampling algorithms including MiDaS by

```
make clean
make
./run_sampling.sh
```

### Sampling

* src/main.cpp
after `make`, type `./bin/Sampling arguments`

```
required arguments:
  --algorithm ALGORITHM       --algo_opt  ALGORITHMOPTION: choose one of the below
              es                          [ global_deg_min / global_deg_max / global_deg_avg ]: for MiDaS_Basic or Random Hyperedge Sampling(RHS) / MiDaS_Basic_Max / MiDaS_Basic_Avg
              ns                          [ global_deg ]: for Random Node Sampling(RNS) and Random Degree Node(RDN)
              rw                          [ rw_c ]: for RandomWalk(RW)
              ff                          [ ff_c ]: for ForestFire(FF)
              mgs                         [ add / exchange / remove ]: for Metropolis Graph Sampling Add / Replace / Delete
              answer                      [ - ]: compute properties of dataset (Degree, Size, Pair Degree, Int. Size, SV, CC, GCC, Density, Overlapness, Diameter)
              find_skewness               [ - ]: compute skewness of the degree distribution of dataset
              test_dynamic                [ - ]: test code for dynamic updating the distributions of node degrees, hyperedge sizes, node-pair degrees, hyperedge-pair intersection sizes
              adjust_time                 [ - ]: calculate the additional elapsed time for computing the D-Statistics in degree distributions
              cal_phi_dist                [ - ]: observe theorem 2
  --dataname DATA: dataset file name
  --inputpath INPUTPATH: path for the directory of the dataset
  --target_portion PORTION: sampling portion
  --eval_opt EVAL: [degree/avg] When running mgs, choose degree if you aim to preserve node degrees. On the other hand, choose avg if you aim to 
                                preserve node degrees, hyperedge sizes, node-pair degrees and hyperedge intersection sizes at the same time.
  --repeat REPEAT: the order of repeating
```
```
optional arguments depending on ALGORITHM:
  ALGORITHM
  es or ns      --alpha ALPHA: the degree of bias towards those with high-degree nodes
  mgs           --turn TURN: the number of iterations for MGS-DT
  rw            --maxlength MAXLENGTH --restart RESTART: MAXLENGTH for the maximum number of steps and RESTART for restart probability
  ff            --p P --q Q: p and q for the parameters same in HyperFF
  test_dynamic  --num_tries NUMTRIES: the number of iterations
```

* MiDaS

```
cd analyze
python midas.py --data DATA --portion SAMPLINGPORTION
```

* Ablation Study (MiDaS-Basic-Max, MiDas-Basic-Avg, MiDas-Basic-NS)

```
cd analyze
python ablation_study.py --select 0 --algotype es --opt global_deg_max    # MiDaS-Basic-Max
python ablation_study.py --select 0 --algotype es --opt global_deg_avg    # MiDaS-Basic-Avg
python ablation_study.py --select 0 --algotype es --opt global_deg_min    # MiDaS-Basic-Grid in the searchspace for ablation study
python ablation_study.py --select 0 --algotype ns --opt global_deg        # MiDaS-Basic-NS
```

### Find Properties

After running `main.cpp`, properties of the sampled hypergraphs are already saved in the directory except for overlapness, singular values and diameter. Thus, we find remaining properties of the sampled hypergraphs and then calculate distances from the entire hypergraph with respect to ten properties.

* analyze_sv.py

```
python analyze.py
default setting:
  find singular values in all results from seventeen sampling methods(including ablation study) in multiple settings(eleven datasets and five sampling portions)
```
```
arguments:
  --repeat_str REPEATLIST: pass REPEATLIST where indexes are separated by "," if you want to consider *only* results from these repeat indexes. (eg. "1,2,3")
  --portion_str PORTIONLIST: pass PORTIONLIST where sampling portions are separated by "," if you want to consider *only* these portions. (eg. "0.30,0.20")
  --dataname DATALIST: pass DATALIST where datasets are separated by "," if you want to consider *only* sampled hypergraphs from these datasets. (eg. "email-Eu-full,coauth-MAG-Geology-full")
  --algoname ALGORITHMLIST: pass ALGORITHMLIST where algorithms are separated by "," if you want to consider *only* results from these algorithms. (eg. "midas,tihs")
```

For large datasets (threads and coauth domain), run `preprocess_sv.py` and use matlab to find singularvalues fast(refer to `script_for_sv_geology.m`)

* helper.py

```
python helper.py --overlapness --diameter # find overlapness and diameter of the sampled hypergraphs
python helper.py --get_eval # calculate D-Statistics and relative difference
python helper.py --repeat_agg # average of repetition
```
```
default setting:
  run over all results from seventeen sampling methods(including ablation study) in multiple settings(eleven datasets and five sampling portions)
```
```
arguments:
  --inputpath DIRPATH: use if you want to consider *only* results in this directory
  --dataname DATANAME: use if you want to consider *only* sampled hypergraphs from this dataset
  --algorithmlist ALGORITHMLIST: pass ALGORITHMLIST where algorithms are separated by "," if you want to consider *only* results from these algorithms. (eg. "midas,tihs")
  --portion PORTIONLIST: pass PORTIONLIST where sampling portions are separated by "," if you want to consider *only* results from these portions. (eg. "0.30,0.20")
```

### Plot Figures

In the *analyze* directory,

* analyze_result.py (Make Table)
* draw_figures.py
* observation.py
* time_eval.py
* ablation_study.py
* theorem_plot.py

You can see how to run these codes in `analyze/run.sh`
- - -

## Environment

The environment of running codes is specified in `environment.yml`
