# Energy Efficiency Maximization for C-RANs: Discrete Monotonic Optimization, Penalty, and $\ell_{0}$-Approximation Methods

This repo is cloned from [this CodeOcean capsule](https://codeocean.com/capsule/4590968/tree/v1), which contains the code for the following scientific paper:

Kien-Giang Nguyen, Quang-Doanh Vu, Markku Juntti, and Le-Nam Tran, "Energy Efficiency Maximization for C-RANs: Discrete Monotonic Optimization, Penalty, and $\ell_{0}$-Approximation Methods," IEEE Trans. Signal Process., vol. 66, no. 17, pp. 4435-4449, Sep. 2018.

## Instructions
The "**main_optimal.m**" script plots the convergence of the proposed branch-and-bound method described in Algorithm 1 in the paper. It makes use of [MOSEK](https://www.mosek.com/) directly (i.e., not through a parser such as CVX or YALMIP) as a convex solver to solve various conic convex problems in the paper. If you download the code and run it on your machine, then you need to manually install MOSEK yourself.
