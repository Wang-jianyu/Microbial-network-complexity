# Microbial-network-complexity
This repository contains code for the manuscript.
## 1._P_-value calculate
To obtain P values, we performed permutation and bootstrap procedures with 1000 iterations each. We then obtained the measure-specific P value as the probability of the null value (the mean of permutation distribution) under a normal distribution curve fitted to the mean and standard deviation of the bootstrap distribution.
The corresponding R script is **01 p value_calculate.R**.
## 2.Network inference
The co-occurrence network were obtained by calculating Spearman correlation and Jaccard dissimilarity measures.
The corresponding R script is **02 network_infer.R**.
