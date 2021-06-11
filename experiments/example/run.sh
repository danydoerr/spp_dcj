#!/bin/sh

# generate ILP
../../scripts/spp_dcj.py -m Example1.idmap -a 0.5 -b 0.25 SpeciesTree.txt \
    AdjacenciesExample1.txt > Example1.ilp 2> Example1.spp_dcj.log
# run ILP
gurobi_cl ResultFile=Example.sol Example1.ilp > Example1.gurobi.log
# extract adjacencies from solution
../../scripts/sol2adjacencies.py Example1.sol Example1.idmap > \
    PredictedAdjacenciesExample1.txt
# visualize candidate and predicted adjacencies
 ../../scripts/visualize_genomes.py -i PredictedAdjacenciesExample1.txt \
     AdjacenciesExample1.txt > PredictedAdjacenciesExample1.pdf
