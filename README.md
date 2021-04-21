# Microbial-network-complexity
This repository contains code for the manuscript.
## 1._P_-value calculate
To obtain P values, we performed permutation and bootstrap procedures with 1000 iterations each. We then obtained the measure-specific P value as the probability of the null value (the mean of permutation distribution) under a normal distribution curve fitted to the mean and standard deviation of the bootstrap distribution.
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
## 5.Stability
We used network robustness as a measure of the stability of interaction networks to species losses. Following this rationale, we simulated three projected scenarios of species extinction: taxa removal in the order of (ⅰ) least-to-most and (ⅱ) most-to-least abundance, and finally (ⅲ) random species removal.The network robustness estimates were then extended to the sub-networks of each sampling plot.
The corresponding R script is **05 stability.R** and **06 sub stability.R**.
