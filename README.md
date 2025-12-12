# A Reduction-Driven Local Search Algorithm For Generalized Independent Set Problem
The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper A Reduction-Driven Local Search Algorithm For Generalized Independent Set Problem by Yiping Liu, Yi Zhou, Zhenxiang Xu, Mingyu Xiao and Jin-Kao Hao.

# Description
We study the generalized independent set (GIS) problem.
In an undirected graph $G=(V,E)$, where $V$ is the vertex set and $E$ is the edge set, each vertex $v\in V$ is assigned a profit $w(v)\in \mathbb{R}$ and each edge $e\in E$ is either a permanent edge or a removable edge. 
If $e$ is a removable edge, then it is further associated with a penalty $p(e)\in \mathbb{R}$. 
The GIS problem asks for a subset $I\subseteq V$ that contains no permanent edges with both endpoints in $I$, and maximizes the difference between the sum of vertex profits in $I$ and the sum of penalties of removable edges with both endpoints in $I$. An example is shown in the following figure. 
The set of red vertices is a solution with a net benefit of $2-1+6+5=12$.

<p align="center">
  <img src="https://github.com/PlutoAiyi/RLS/blob/main/GIS%20example.png?raw=true"  alt="Sublime's custom image" />
</p>


The GIS problem has a number of applications in modern areas such as forest harvesting, competitive facility location, social network analysis, and even machine learning. However, solving the GIS problem in large-scale real-world networks is still computationally intractable. 
To address the issue, we first propose 14 reduction rules that can reduce the input graph with optimality guarantees.
We then present a reduction-driven local search (RLS) algorithm that integrates these
reduction rules into the pre-processing, the initial solution generation, and the local
search operation phases in a computationally efficient way. We further evaluate RLS with 278 graphs arising from different application scenarios. The performance of RLS
is highly competitive – It finds much better solutions than other known solvers for a majority of graphs, and it delivers high-quality solutions for graphs with more than 10
million vertices, while every known approach fails. Analysis also reveals that the data reduction plays a key role in achieving such a competitive performance.

This project contains a dataset folder which contains the datasets and the link to datasets used in the paper.


# Prerequisites
The codes are implemented under Ubuntu 22.04. Boost C++ 17 is also required for running the codes.

# Build and Run
In Linux/Ubuntu operating system, execute the following command to build the program.
```
make
```
To run the executable file, call the following command.

```
./RLS datasets/[dataset name.txt]
```
Be sure to clean all the dependencies and executable files before building a different version of the code.


# Note
The GIS.mod file is used for CPLEX v22.12 to exactly solve the instances.
We also provide the link for benchmark algorithms:

**CB&B**:The exact branch-and-bound algorithm using the UBQP formulation proposed by Hosseinian and Butenko (2019). The source code was implemented in C++ and made available by the authors. We share the source code in https://drive.google.com/file/d/1h7WCpe3mz-WKQftBW8l3EipsIhjpJRPB/view?usp=sharing

**LA-B&B**: The exact branch-and-bound algorithm using a Lagrange relaxation technique by Zheng et al. (2024). The source code is implemented in JAVA and made public in https://github.com/m2-Zheng/GISP.

**ILS-VND**: The iterative local search algorithm proposed by Nogueira et al. (2021). https://sites.google.com/site/nogueirabruno/software

**ALS**: The adaptive multi-neighborhood local search algorithm proposed by Zheng et al. (2024). The source code is implemented in JAVA and made public in https://github.com/m2-Zheng/GISP.

