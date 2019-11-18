# TropicalOT

This repository provides code to implement the numerical experiments in [Lee et al.](https://arxiv.org/abs/1911.05401), which solves the optimal transport problem on the ambient space of phylogenetic trees (i.e., the tropical projective torus).  Specifically, the Wasserstein-1 and 2 distances are calculated, with the tropical metric as ground metric.

Wasserstein distances are metrics on spaces of probability measures.  Intuitively, they measure the effort required to recover the probability mass of one distribution in terms of an efficient reconfiguration of another.  In terms of optimal transport, they yield the solution to the optimal transport problem when the cost function of moving a mass from one location to another is no more than the distance between the locations.  

For more detail on the underlying theory of optimal transport and Wasserstein distances on the tropical projective torus , see [Lee et al.](https://arxiv.org/abs/1911.05401) and the references therein.

This repository consists of two separate directories `Tropical_Wasserstein1` and `Tropical_Wasserstein2` consisting of C++ code to calculate the tropical Wasserstein-1 and 2 distances.

## Running TropicalOT
### Requirements
The following components are required to run the code in this repository:
* g++ is the C++ compiler for GNU compiler collection
* fftw3 is a C subroutine to compute discrete Fourier transform, available at [fftw.org](http://www.fftw.org/)
* python 3 is available at [python.org](https://www.python.org/downloads/) and is used to plot output (optional)

Open a terminal (command line prompt if using Windows) and type the following:
```
g++ -O3 tropical_w1.cpp -o main.exe -lfftw3
```
Replace `tropical_w1.cpp` with `tropical_w2.cpp` as appropriate to compute tropical Wasserstein-2 distances. 

Once the compilation is complete, the following command runs the code for both the tropical Wasserstein-1 and 2 distances (with the appropriate script called in the previous line):
```
	./main.exe [n1] [n2] [h] [tau] [tolerance] [max_iteration]
```
Here, n1 and n2 represent the discretization of the 2D square domain [0,1]x[0,1], h and tau represent step sizes of the G-Prox PDHG algorithm, and tolerance and max_iteration are specified for the iterations of the G-Prox PDHG algorithm.

Here is an example:
```
	./main.exe 64 64 16 0.1 0.5 0.00001 1000
```

Running the command line for the tropical Wasserstein-1 distance will produce the following five `.csv` files in the `Data` directory:
* mx.csv
* my.csv
* rho0.csv
* rho1.csv
* Parameters.csv

These may be plotted by typing the following command in a terminal, which will save the plots in the `Images` directory:
```
	python plot.py
```

Similarly, running the command line for the tropical Wasserstein-2 distance will produce the following two `.csv` files in the `Data` directory:
* rho.csv
* parameters.csv

A video for these images may be created by typing the following command in a terminal, which will create and save a file named `video.mp4`

## Relevant References
*Tropical Optimal Transport and Wasserstein Distances in Phylogenetic Tree Space*
Wonjun Lee, Wuchen Li, Bo Lin, and Anthea Monod
[arxiv:1911.05401](https://arxiv.org/abs/1911.05401)

## Contact
We appreciate any questions or feedback you have on our code.  Please contact [Wonjun Lee](mailto:wlee@math.ucla.edu) or [Anthea Monod](mailto:antheam@tauex.tau.ac.il).
