# Successive continuation based optimization

The necessary conditions for constrained optimization problems include the original system (constraints and objective function(al)) and an adjoint system. Such necessary conditions are typically formulated as a nonlinear boundary-value problem (BVP). With discretization applied, the BVP is converted into a set of nonlinear equations. Constructing reasonable initial guesses for solving such nonlinear equations is challenging. Kernévez and Doedel [1] proposed a succssive continuation scheme to tackle the challenge, in which solutions to the necessary conditions for locally optimal solutions are found at the end of a sequence of easily initialized separate stages of continuation. In particular, the first run is initialized with trivial Lagrange multipliers.

In [2], the authors etasblished staged construction of adjoint system and original system. A predefined library for algebraic, differential and integral operators and their adjoints was developed. This library enables automated generation of the necessary conditions. Such conditions are satisfied using the sucessive continuation method. This functionality was implemented with the November
2017 release of COCO, a MATLAB-based open-source toolbox for numerical continuation. Please refer 
https://sourceforge.net/projects/cocotools/ for more info of COCO.

In [3], the authors generalized the successive continuation scheme proposed in [1] to the case of simultaneous equality and inequality constraints. A key enabler of the proposed generalization is the use of complementarity functions to define relaxed complementary conditions, followed by the use of continuation to arrive at the limit required by the KKT theory. This functionality is avaliable in the March 2020 release of COCO.

This repository presents the code for the examples in [2] and [3].

Read More: https://epubs.siam.org/doi/abs/10.1137/17M1143563; https://doi.org/10.1016/j.amc.2020.125058

## References
[1] Kernévez, J. P., & Doedel, E. J. (1987). Optimization in bifurcation problems using a continuation method. In Bifurcation: Analysis, Algorithms, Applications (pp. 153-160). Birkhäuser Basel.

[2] Li, M., & Dankowicz, H. (2018). Staged construction of adjoints for constrained optimization of integro-differential boundary-value problems. SIAM Journal on Applied Dynamical Systems, 17(2), 1117-1151.

[3] Li, M., & Dankowicz, H. (2020). Optimization with equality and inequality constraints using parameter continuation. Applied Mathematics and Computation, 375(15), art.no. 125058.
