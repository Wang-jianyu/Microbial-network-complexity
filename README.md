# Microbial-network-complexity
This repository contains code for the manuscript.
## 1._P_-value calculate
To obtain _P_ values, we performed permutation and bootstrap procedures with 1000 iterations each. We then obtained the measure-specific _P_ value as the probability of the null value (the mean of permutation distribution) under a normal distribution curve fitted to the mean and standard deviation of the bootstrap distribution.
The corresponding R script is **01 p value_calculate.R**.
## 2.Network inference
The co-occurrence network were obtained by calculating Spearman correlation and Jaccard dissimilarity measures.
The corresponding R script is **02 network_infer.R**.
## 3.Network property
We calculate main topological properties potentially relevant for community roles and functioning.
The corresponding R script is **03 network_property.R**.
## 4.Subnetwork inference
The obtained site-level network meta-matrix were then used to sub-set networks matrices for each sampling plot by preserving OTUs presented within the plot and all edges among them in the site-level network.
The corresponding R script is **04 sub_net_infer.R**.
