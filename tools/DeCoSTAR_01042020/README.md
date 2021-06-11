# DeCoSTAR

DeCoSTAR is a program designed to study the evolution of links between genes (i.e., adjacencies) [0].

Given a species tree S, a set of gene family (unrooted) tree distributions G, a set of extant adjacencies A and a set of costs for adjacencies events (adjacency gain and adjacency breakage), DeCoSTAR can compute:

    (1) reconciled gene trees R from the gene trees in G, such that these forms a most parsimonious reconciliation between S and G according to the TERA algorithm described in [1].
    (2) an history of the given adjacencies in A along the reconciled gene trees R such that this history minimizes an adjacency gains and breakages cost(with respect to their relative costs); according to the models described in [2][3].

The history of adjacencies comes in the form of one or several adjacency trees (phylogenetic trees in which nodes represents adjacencies between the nodes of gene trees). We refer to such an history as an adjacency forest. Note that DeCoSTAR can, but usually does not output adjacency forest but instead provide lists of inferred adjacencies in ancestral species.

Rather than an adjacency forest minimizing the number of adjacency gains and breakages, DeCoSTAR can sample adjacency forests in such a manner that adjacency forests with a lower adjacency cost have a higher chance to be sampled. This is done according to the algorithm described in [4], extended to include transfers.

It is also possible to use DeCoSTAR to infer adjacencies in extant species as described in [5] by using the scaffolding mode.

Installation
============

DeCoSTAR requires:
 * Bio++ ( http://biopp.univ-montp2.fr/) 
 	* installation instructions : http://biopp.univ-montp2.fr/wiki/index.php/Installation#Getting_the_source_files
 	* the current version has been succesfully build using the following commits:
		* bpp-core : df2daf3b70ebd74d1b4d9c318d8db938228a349a
		* bpp-phyl : 630247204cc7c0d85cf2cb90af3a1a5feb0dac86
		* bpp-seq  : 631d58e919286bfb40d05bc63c201b31405f6d35

 * boost (http://www.boost.org/)


To install DeCoSTAR, place yourself in the root directory of this repository.

First, edit the makefile lines:

```
BPP_INCLUDE=$(HOME)/local/bpp/dev/include
BPP_LIB=$(HOME)/local/bpp/dev/lib

BOOST_INCLUDE=/usr/include
BOOST_LIB=/usr/lib
```

to make them reflect where you have installed bio++ and boost

type :
```
make bin/DeCoSTAR
```


References
==========

[0] Duchemin W., Anselmetti Y., Patterson M., Ponty Y., Berard S., Chauve C., Scornavacca C., Daubin V., Tannier E. (2017): DeCoSTAR: Reconstructing the Ancestral Organization of Genes or Genomes Using Reconciled Phylogenies. In: Genome biology and evolution, vol. 9 pp.1312-1319.

[1] Celine Scornavacca, Edwin Jacox, and Gergely Szöllősi. Joint Amalgamation of Most Parsimonious Reconciled Gene Trees. Bioinformatics (2014): btu728.

[2] Sèverine Bérard, Coralie Gallien, Bastien Boussau, Gergely J. Szöllősi, Vincent Daubin and Eric Tannier. Evolution of gene neighborhoods within reconciled phylogenies. Bioinformatics (Oxford, England) Vol. 28, No. 18 (2012) p. i382-i388

[3] Murray Patterson, Gergely Szöllősi, Vincent Daubin, Eric Tannier. Lateral gene transfer, rearrangement, reconciliation. BMC Bioinformatics Vol. 14, Suppl. 15 (2013) S4

[4] Cedric Chauve, Yann Ponty and João Paulo Pereira Zanetti. Evolution of genes neighborhood within reconciled phylogenies: an ensemble approach. BMC Bioinformatics Vol. 16 Suppl. 19  (2015) S6

[5] Yoann Anselmetti, Vincent Berry, Cédric Chauve, Annie Chateau, Eric Tannier and Sèverine Bérard. Ancestral gene synteny reconstruction improves extant species scaffolding. BMC genomics Vol. 16, Suppl. 10 (2015) S11  
