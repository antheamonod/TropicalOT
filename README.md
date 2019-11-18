# TropicalOT

This repository provides code to implement the numerical experiments in [Lee et al.](https://arxiv.org/abs/1911.05401), which solves the optimal transport problem on the ambient space of phylogenetic trees (i.e., the tropical projective torus).  Specifically, the Wasserstein-1 and 2 distances are calculated, with the tropical metric as ground metric.

Wasserstein distances are metrics on spaces of probability measures.  Intuitively, they measure the effort required to recover the probability mass of one distribution in terms of an efficient reconfiguration of another.  In terms of optimal transport, they yield the solution to the optimal transport problem when the cost function of moving a mass from one location to another is no more than the distance between the locations.  For more detail on the underlying theory, see [Lee et al.](https://arxiv.org/abs/1911.05401) and the references therein.

This repository consists of two separate directories `Tropical_Wasserstein1` and `Tropical_Wasserstein2` consisting of C++ code to calculate the tropical Wasserstein-1 and 2 distances.
