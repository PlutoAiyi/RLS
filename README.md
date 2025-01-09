# A Reduction-Driven Local Search Algorithm For Generalized Independent Set Problem
The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper A Reduction-Driven Local Search Algorithm For Generalized Independent Set Problem by Yiping Liu, Yi Zhou, Zhenxiang Xu, Mingyu Xiao and Jin-Kao Hao.

# Description
We study the \textit{generalized independent set} (GIS) problem.
In an undirected graph $G=(V,E)$, where $V$ is the vertex set and $E$ is the edge set, each vertex $v\in V$ is assigned a profit $w(v)\in \mathbb{R}$ and each edge $e\in E$ is either a \textit{permanent edge} or a \textit{removable edge}. 
If $e$ is a removable edge, then it is further associated with a penalty $p(e)\in \mathbb{R}$. 
The GIS problem asks for a subset $I\subseteq V$ that contains no permanent edges with both endpoints in $I$, and maximizes the difference between the sum of vertex profits in $I$ and the sum of penalties of removable edges with both endpoints in $I$. An example is shown in the following figure.
![alt text](https://github.com/PlutoAiyi/RLS/blob/main/GIS%20example.png?raw=true)

The set of red vertices is a solution with a net benefit of $2-1+6+5=12$.



has a number of applications in modern
areas such as forest harvesting, competitive facility location, social network analysis,
and even machine learning. However, solving the GIS problem in large-scale real-
world networks is still computationally intractable. To address the issue, we first
propose 14 reduction rules that can reduce the input graph with optimality guarantees.
We then present a reduction-driven local search (RLS) algorithm that integrates these
reduction rules into the pre-processing, the initial solution generation, and the local
search operation phases in a computationally efficient way


The second-best congestion pricing (SBCP) problem is one of the most challenging problems in transportation due to its two-level hierarchical structure. In spite of various intriguing attempts for solving SBCP, existing solution methods are either heuristic without convergence guarantee or suitable for solving SBCP on small networks only. In this paper, we first reveal some convexity-based structural properties of the marginal value function reformation of SBCP and then, by effectively exploiting these structural properties, we propose two dedicated decomposition methods for solving SBCP on large-scale networks which are different from existing methods in that they avoid linearizing nonconvex functions. We establish the convergence of the two decomposition methods under commonly used conditions and provide the maximum number of iterations for deriving an approximate stationary solution. The computational experiments based on a collection of real road networks show that in comparison with three existing popular methods, the two proposed methods are capable of solving SBCP on larger-scale networks; and for instances that can be solved by existing methods, the two proposed methods are substantially faster.

This project contains two folders: PD, RPD.

PD: This folder contains the data, source codes, makefile and results of the penalty decomposition methods.
RPD: This folder contains the data, source codes, makefile and results of the relaxed penalty decomposition methods.
