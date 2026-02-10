#  Outline - Networks as a Proxy for Allosteric Signaling in PDZ

## Introduction

- We know of allostery's existence, but don't know how to measure it
- This paper aims to provide a method to understand and measure allostery using network theory.

## Background

### Networks

- Networks have been used to model proteins before
- $G = (V, E)$ is a graph, where $V$ is the set of vertices and $E$ is the set of edges
- $G$ is a weighted directed graph, where each edge has a weight
- $G$ is a sparse matrix, where each element is the weight of the edge between two vertices
- $G$ is a symmetric matrix, where each element is the weight of the edge between two vertices
- $G$ is a directed graph, where each edge has a weight
- $G$ is a weighted directed graph, where each edge has a weight
- $G$ is a sparse matrix, where each element is the weight of the edge between two vertices
- $G$ is a symmetric matrix, where each element is the weight of the edge between two vertices
- $G$ is a directed graph, where each edge has a weight

### Random Walks

### Markov

- Definition of a Markov Chain
  - A Markov chain is an absorbing Markov Chain if:
    - It has at least one absorbing state and it is possible to reach an absorbing state from any state
- Absorption Probabilities
    - Start from a symmetric correlation/weight matrix and take absolute values
    - Set the diagonal entries to zero (no self-transitions), $p_{ii} = 0$ for all $i$
    - Apply a high-correlation threshold (0.85): values below the threshold are treated as extremely small (epsilon) edges to ensure all transient states can reach an absorbing state
    - Normalize each row so that outgoing probabilities sum to 1, yielding a Markov transition matrix
    - Make target nodes absorbing by zeroing their outgoing probabilities and setting their self-loop probability to 1, $p_{ii} = 1$ for all $i \in$ absorbing nodes
    - Partition the transition matrix into $Q$ (transient → transient) and $R$ (transient → absorbing), compute the fundamental matrix $F = (I − Q)^{-1}$, and obtain the absorption probabilities as $B = F R$, where $B_{ij}$ is the probability that a walk starting at transient state $i$ is eventually absorbed in absorbing state $j$
- TODO: Get absorption probabilities with absorbing nodes as A-site and L-site

### Bio stuff

## Methods

### Trajectory Data


### Software packages

- Bio
- biopython
- matplotlib
- mdtraj
- networkx
- numpy
- pandas
- plotly
- scipy
- seaborn
- ydata_profiling

## Results

### My Method

### Setup system

- Trajectories -> DCCM -> Network -> Random Walks -> Absorption Probabilities -> Bottleneck Analysis

#### DCCM

- Dynamic Cross-Correlation Matrix
- Computed using AMBER's cpptraj tool
- ```
  parm $prmtop_file
  trajin $mdcrd_file
  reference $mdcrd_file 1
  
  box auto
  center @CA
  rms reference mass out rmsd_${base_name}.dat
  matrix correl @* out $output_file
  
  run
  EOF
  ```
    - TODO: have a sentence explaining why basing off of frame 1. Check Lakhani's paper.
    - rms reference... aligns trajectory to reference for all atoms, removing global rotations and translations from the correlation matrix
$$
M_{ij} =
\frac{\langle \mathbf{V}_i \cdot \mathbf{V}_j \rangle -
\langle \mathbf{V}_i \rangle \cdot \langle \mathbf{V}_j \rangle}
{\sqrt{
(\langle \mathbf{V}_i^2 \rangle - \langle \mathbf{V}_i \rangle^2) \
(\langle \mathbf{V}_j^2 \rangle - \langle \mathbf{V}_j \rangle^2)
}}
$$

#### Network Creation

- DCCM is an adjacency matrix, where each element is the correlation between two atoms
- Threshold of |0.8| to simplify network and only look at high correlations edges
- $e_{ij} = \begin{cases} |M_{ij}| & \text{if } |M_{ij}| > 0.8 \\ 0 & \text{otherwise} \end{cases}$
- 
- TODO: redo random walks completely, with the different absorption nodes (A-site and L-site). Use the transition matrix matricies, not the transiion absorption matrices.

### Random Walk Analyses

- ```python
  def random_walk(P: sp.csr_matrix, start: int, max_steps: int, target_nodes: set[int], rng) -> tuple[list[int], Literal[1, 2, 3, 4], list[float]]:
    ```
    - P is the transition matrix
    - start is the starting node
    - max_steps is the maximum number of steps to take in the walk
    - target_nodes is the set of target nodes to terminate the walk at
    - rng is the random number generator
    - Returns a tuple containing the path taken during the walk, the termination code, and the list of log-transformed probabilities of each step in the walk

### Absorption Probabilities

- Threshold of |0.85|

### Bottleneck Analysis

## Discussion

## Conclusion

## References