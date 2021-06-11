# SPP_DCJ

## Introduction
The Small Parsimony Problem (SPP) aims at finding the gene orders at internal nodes of a given phylogenetic tree such that the overall genome rearrangement distance along the tree branches is minimized. This problem is intractable in most genome rearrangement models, especially when gene duplication and loss are considered.
`SPP_DCJ` is _Integer Linear Program_-based algorithm to solve the SPP for natural genomes, _i.e._, genomes that contain conserved, unique, and duplicated markers. The evolutionary model that we consider is the _DCJ-indel model_ that includes the Double-Cut and Join rearrangement operation and the insertion and deletion of genome segments. 

`SPP_DCJ` is an extension of [DING](https://gitlab.ub.uni-bielefeld.de/gi/ding).

## Dependencies
- [Gurobi](https://www.gurobi.com/products/gurobi-optimizer/) or [CPLEX](https://www.ibm.com/analytics/cplex-optimizer)
- [Python 3](https://www.python.org/downloads/)
- Python 3 libraries:
    - Matplotlib
    - NetworkX
    - Numpy
    - Pandas
    - Ete3

## Optional Dependencies
- Python 3 libraries
    - Snakemake

## How to run 

The following steps show howto run `SPP_DCJ` with Gurobi.

0. `SPP_DCJ` requires (i) a given phylogeny and (ii) a table with candidate adjacencies for all genomes corresponding to nodes of the given phylogeny. The candidate adjacency table has the following columns:

    ```#Species        Gene_1  Ext_1   Species Gene_2  Ext_2   Weight```

1. Generate ILP (`-a` and `-b` are parameters of the objective function, set here to 0.5 and 0.25, respectively):
    
    ```scripts/spp_dcj.py -m experiments/example/Example1.idmap -a 0.5 -b 0.25 experiments/example/SpeciesTree.txt experiments/example/AdjacenciesExample1.txt > experiments/example/Example1.ilp 2> experiments/example/Example1.spp_dcj.log```

2. Run ILP
    
    ```gurobi_cl ResultFile=experiments/example/Example.sol experiments/example/Example1.ilp > experiments/example/Example1.gurobi.log```

3. Extract adjacencies from solution

    ```scripts/sol2adjacencies.py experiments/example/Example1.sol experiments/example/Example1.idmap > experiments/example/PredictedAdjacenciesExample1.txt```

4. Visualize candidate and predicted adjacencies

    ```scripts/visualize_genomes.py -i experiments/example/PredictedAdjacenciesExample1.txt experiments/example/AdjacenciesExample1.txt > experiments/example/PredictedAdjacenciesExample1.pdf```

## Example
An the code above in included in the example in `experiments/example`.



