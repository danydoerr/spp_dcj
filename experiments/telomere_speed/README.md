# HOWTO run the telomere_speed workflow

1. edit the config and adjust parameters:
    - path to GUROBI binary `gurobi_bin: /path/to/binary` (unless it is available in the system
      path)
    - path to DING scripts directory `ding_dir:        python2 /path/to/ding/scripts`. DING is available [here](https://gitlab.ub.uni-bielefeld.de/gi/ding).
2. *(optional)* adjust simulation parameters:
    - all variables maintain *lists* of of parameter choices, so that multiple
      choices can be run simultaneously
3. call `snakemake`
