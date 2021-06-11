# HOWTO run the simulation workflow

1. edit the config and adjust parameters:
    - path to GUROBI binary `gurobi_bin: /path/to/binary` (unless it is available in the system
      path)
2. *(optional)* adjust simulation parameters:
    - all variables maintain *lists* of of parameter choices, so that multiple
      choices can be run simultaneously
3. call `snakemake`
