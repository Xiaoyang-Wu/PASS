# PASS
Codes for reproducing the numerical results in *Optimal Subsampling via Predictive Inference (2023)* by Xiaoyang Wu, Yuyang Huo, Haojie Ren and Changliang Zou.

# Description of PASS
In the big data era, subsampling or sub-data selection techniques are often adopted to extract a fraction of informative individuals from the massive data. Existing subsampling algorithms focus mainly on obtaining a representative subset to achieve the best estimation accuracy under a given class of models. In this paper, we consider a semi-supervised setting wherein a small or moderate sized “labeled” data is available in addition to a much larger sized “unlabeled” data. The goal is to sample from the unlabeled data with a given budget to obtain informative individuals that are characterized by their unobserved responses. We propose an optimal subsampling procedure
that is able to maximize the diversity of the selected subsample and control the false selection rate (FSR) simultaneously, allowing us to explore reliable information as
much as possible. The key ingredients of our method are the use of predictive inference for quantifying the uncertainty of response predictions and a reformulation of the objective into a constrained optimization problem.

# Contents

- **R-file**: Algorithm classes and function definitions.
- **Simulation Codes**: Codes for reproducing the illustration example plot and the simulated results.
- **Real-data**: Codes for reproducing the real-data experiment and the data set used.
