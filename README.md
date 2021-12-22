# Private Optimization

Code to accompany, and generate figures for, the [working paper](https://arxiv.org/abs/2103.11003) "Differentially private inference via noisy optimization" M. Avella-Medina, C. Bradshaw and P.L. Loh (2021)

**figs** contains code for generating figures and tables from the paper

**src** contains implementations of noisy gradient descent and noisy newton, for linear and logistic regression. Functions in this main folder satisfy mu-GDP when "private=T" is specified.

**src/extra** contains optimization functions modified to produce certain figures from the paper. Some of those modifications cause the functions to violate the differential privacy guarantee when "private=T" is specified. This is for illustrative purposes only. Code in this folder is not appropriate for real differentially private data analysis tasks.
