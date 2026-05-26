![](schematic.png)

# Overview
This repository contains code and data for generating the primary results from "Hierarchical Community Structure of the Adult Drosophila Connectome Reveals Conserved Circuit Archetypes."

## To run analyses
1. Download and unzip this repository.
2. Navigate to this [zenodo repository]{https://zenodo.org/records/20399924} and download the connectome data and the hierarchical community partition. Place these files in the `data/` directory.
3. Open MATLAB, navigate to the `m/` directory and execute `run_analyses.m`

## Pre-computed derivatives
Note: that `m/run_analyses.m` performs analyses to generate Figures 2-8 (note that Figure 1 is simply a schematic). The scripts for each individual figure (e.g. `figure*.m`) include some computationally intense steps. For convenience, we have run these scripts and saved intermediate outputs so that the most computationaly intense steps are unnecessary. These are available in the `outputs/` directory. However, if you would like to run the scripts in their entirety, you can edit them so that, instead of loading in the pre-generated outputs, you compute them yourself.

We also include code (`py/sample_posterior_distribution_no_threshold.ipynb`) for estimating communities based on the [graph-tool]{https://graph-tool.skewed.de/} library. This procedure is computationally intense and will take a long time, owing to the size of the connectome. If the goal is to simply reproduce the results of the main text, we strongly suggest running `run_analyses.m` using the partition from the [Zenodo dataset]{https://zenodo.org/records/20399924}.

## References
If you use data from this paper, please abide by FlyWire citation guidelines and cite the relevant [primary resources]{https://codex.flywire.ai/about_flywire}. Additionally, please cite our paper:

Betzel, R., Del Rio, O., Labora, N., Dvali, S., Larsen, B., Lynn, C. W., ... & Seguin, C. (2025). Hierarchical Community Structure of the Adult Drosophila Connectome Reveals Conserved Circuit Archetypes. bioRxiv, 2025-12.

All code and data are provided as is. If you have specific questions, please reach out to Rick Betzel (rbetzel@umn.edu).
