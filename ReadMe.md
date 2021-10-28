# MiDaS: Representative Sampling from Real-world Hypergraphs

<img src="correction/correction.jpg" width="900"/>

- - -

We provide source code for, 

(1) **Sampling Hypergraphs**
* *Simple and Intuitive Approaches* : six intuitive approaches(RNS, RDN, RW, FF, RHS, TIHS) experimenting on 11 datasets with a focus on their limitations in preserving 10 properties
* *MiDaS-Basic* : a preliminary version of MiDaS which can improve the limitations of intuitive approaches by proper alpha setting
* *MiDaS*: a full-fledged version which can automatically tune alpha

(2) **Evaluating the Quality of Sub-Hypergraph**

Measure how precisely the sub-hypergraph preserves the structural properties of the entire hypergraph with respect to,

* *node-level statistics*: the distributions of node degrees and node-pair degrees
* *hyperedge-level statistics* : the distributions of hyperedge sizes and intersection sizes
* *graph-level statistics* : average clustering coefficient, density, overlapness, and effective diameter


## Datasets

In the papaer, we used datasets after removing duplicated hyperedges. We preprocessed thirteen datasets collected We preprocessed datasets collected by [Austin R. Benson][https://www.cs.cornell.edu/~arb/]. The datasets used in the paper are available in the "dataset" folder.

## How to Run

#### Example
You can run all thirteen sampling algorithms including MiDaS by

```
./run_sampling.sh
```

#### src/main.cpp


#### analyze


## Environment

The environment of running codes is specified in `environment.yml`
