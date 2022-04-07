# C - Causal Inference

This document describes how to run module B (Causal Inference) of the RAPToR pipeline.

Sequentially run scripts 1-4:

* <b>1-createBlacklist.R</b> - Restrict edges as described on the figure at the bottom.
<br> The algorithm takes as input a matrix "PTR_species_filtered_metadata_major_AMR.csv" 
where each raw is a sample and each column is a random variable;
* <b>2-CausalGraph</b>  - apply constrained-based causal structural learning algorithm (PC-Stable);
* <b>3-writeStyle</b> - apply syling rules to the learned causal structure and transform it to the xgmml format;
* <b>4-ComputeCorrelations</b> - compute Spearman correlation between variables connected with an arc.

<img src="https://github.com/stebliankin/RAPToR/blob/master/images/restrictions.png" width=400>

### Output
* <i>blacklist.csv</i> - list of restricted edges;
* <i>bn_correlation.csv</i> - strength of Spearman correlation for each of the edge;
* <i>boot_strength.csv</i> - bootstrap score for each of the edge;
* <i>directed.csv</i> - all inferenced directed edges;
* <i>undirected.csv</i> - all infereced undirected edges.
