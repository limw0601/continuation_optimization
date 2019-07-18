# Sucessive continuation based optimization

The necessary conditions for constrained optimization problems include original system (constraints and objective function(al)) and adjoint system. Such necessary conditions are typically formulated as nonlinear boundary-value problems (BVP). After discretization, the BVP is converted into a set of nonlinear equations. Constructing reasonable initial guesses for such nonlinear equations is challenging. Kernévez and Doedel [1] proposed a succssive continuation scheme to tackle the challenge, in which solutions to the necessary conditions for locally optimal solutions are found at the end of a sequence of easily initialized separate stages of continuation. In particular, the first run is initialized with trivial Lagrange multipliers.

In [2], the authors etasblished staged construction of adjoint system in parallel with original system. A predefined library for algebraic, differential and integral constraints was developed. The automatically generated necessary conditions are then solved using the sucessive continuation method. The toolbox is avaliable in 2017 Nov release of COCO, a matlab-based toolbox for numerical continuation. Please refer 
https://sourceforge.net/projects/cocotools/ for the info and installation of COCO.

This repository presents the code for the examples in [2]. Code for other interesting optimization problem is included as well.

Read More: https://epubs.siam.org/doi/abs/10.1137/17M1143563

## Reference
[1] Kernévez, J. P., & Doedel, E. J. (1987). Optimization in bifurcation problems using a continuation method. In Bifurcation: Analysis, Algorithms, Applications (pp. 153-160). Birkhäuser Basel.
[2] Li, M., & Dankowicz, H. (2018). Staged Construction of Adjoints for Constrained Optimization of Integro-Differential Boundary-Value Problems. SIAM Journal on Applied Dynamical Systems, 17(2), 1117-1151.
