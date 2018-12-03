# Demo---Deterministic-consensus-maximization-with-biconvex-programming

An efficient and deterministic consensus maximization method called IBCO (as the abbreviation of Iterative Bi-Convex Optimization). 

Published in ECCV 2018 as **oral** presentations.

![alt text](https://github.com/ZhipengCai/ZhipengCai.github.io/blob/master/papers/ECCV18_IBCO.png " ")

About
=====

Consensus maximization is an effective tool for robust fitting in computer vision. This repository contains the demo for IBCO and other prevalent consensus maximization methods (see below for the list). 

If you want to try some methods for your own problem, some personal suggestions are provided at the end of the method list :)

Please refer to the [paper](https://arxiv.org/pdf/1807.09436.pdf) for more details.

The demo is free for non-commercial academic use. Any commercial use is strictly 
prohibited without the authors' consent. Please acknowledge the authors by citing:

```
@inproceedings{cai2018deterministic,
  title={Deterministic Consensus Maximization with Biconvex Programming},
  author={Cai, Zhipeng and Chin, Tat-Jun and Le, Huu and Suter, David},
  booktitle={European Conference on Computer Vision},
  pages={699--714},
  year={2018},
  organization={Springer}
}
```
in any academic publications that have made use of this package or part of it.

------------------------
Contact
------------------------

Homepage:https://zhipengcai.github.io/

Email: zhipeng.cai@adelaide.edu.au

Do not hesitate to contact the authors if you have any question or find any bugs :)

------------------------
List of included methods
------------------------

**Random methods**:

1. [RANSAC](http://delivery.acm.org/10.1145/360000/358692/p381-fischler.pdf?ip=129.127.229.14&id=358692&acc=ACTIVE%20SERVICE&key=65D80644F295BC0D%2E001A23AA3BABC648%2E4D4702B0C3E38B35%2E4D4702B0C3E38B35&__acm__=1543556593_784052ca099a175d04afeade036d626c)

2. [PROSAC](https://ieeexplore.ieee.org/document/1467271#full-text-section)

3. [Guided MLESAC](http://www.robots.ox.ac.uk/~lav/Papers/tordoff_murray_tpami2005/tordoff_murray_tpami2005.pdf)

4. [Locally Optimized RANSAC](http://cmp.felk.cvut.cz/~matas/papers/chum-dagm03.pdf) (LO-RANSAC)

5. [Fixing LO-RANSAC](http://cmp.felk.cvut.cz/software/LO-RANSAC/Lebeda-2012-Fixing_LORANSAC-BMVC.pdf)

6. [USAC](http://people.inf.ethz.ch/pomarc/pubs/RaguramPAMI13.pdf)

**Deterministic methods**:

7. [The Exact Penalty (EP) method](https://arxiv.org/pdf/1710.10003.pdf)

8. [The Smooth Surrogate (SS) method](https://link.springer.com/content/pdf/10.1007/978-3-319-78199-0_21.pdf)

9. [The ADMM-based method](http://bmvc2018.org/contents/papers/0568.pdf) (to be added)

10. [IBCO](https://arxiv.org/pdf/1807.09436.pdf) 

**(personal) Tips for method choosing**: 

+ The runtime of random methods is exponential to 1) the proportion of outliers and 2) the dimension of the problem. If your problem does not have much outliers (say < 50%) and a large dimension (say within 10), random methods are very efficient and usually provide reasonable results. Try method 4~6 in this case.

+ The runtime of deterministic methods is polynomial to the size of input data. If your data size is not extremely large, deterministic methods are fast even with high outlier rates. Try SS (method 8) if you want the highest efficiency. IBCO (method 10) usually returns the best solution among all deterministic methods.

+ Using random methods to initialize deterministic methods is an effective strategy. Deterministic methods are sensitive to initializations, but usually perform better than random methods with a reasonable initial solution. Random methods are capable of returning a reasonable solution, regardless of the initialization. Try LO-RANSAC + IBCO. IBCO is highly suitable for refining the solution of other local methods. Usually it is capable of increasing the consensus size of RANSAC variants by more than 10%. 

+ In practice, if you want high efficiency for problems with high outlier rates or high dimensions, you can terminate random methods earlier and give the current best solution to IBCO. 

-----------------------------------------
List of addressed problems in the demo
-----------------------------------------

**Linear**:

1. Linear regression

**Nonlinear**:

2. Homography estimation

3. Triangulation (to be added)

**With non-convex constraints**:

4. Fundamental matrix estimation (the rank-2 constraint)

IBCO - Introduction
===================

IBCO bisects over all possible consensus sizes and returns the largest one such that the decision version of the consensus maximization problem (briefed as "decision problem" later on) is feasible. 

The major contribution of IBCO is a "lifted" decision problem, which is equivalent to the original problem and bi-convex. Using this lifted form, we can solve each decision problem efficiently by standard bi-convex optimization.

Getting Started
===============

The demo has been tested on Linux (Ubuntu 14.04 LTS 64-bit). Except USAC whose source code is in C++, all other methods are implemented using MATLAB 2017a. 

Due to different languages, we currently do not include USAC in the main demo, while the source code is provided.

-------------
Run the demo
-------------

1. Clone this repository. 

2. Run function "demo()" in MATLAB.

Please refer to "demo.m" file for detailed code explanations.

------------------------
(Optional) Enable Gurobi
-------------------------

This demo uses 'sedumi' as the default linear programming solver, which is used in EP and SS. [Gurobi](http://www.gurobi.com/) was used instead in our experiments, which is faster than 'sedumi'. 

If you have a valid license, you can follow the instruction in the "demo.m" file to switch 'sedumi' to 'Gurobi'.


