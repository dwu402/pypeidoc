---
title: Pypei Variants
author: dwu402
header-includes: |
    \usepackage[section]{placeins}
---

We document here the behaviour of variants of pypei on a testbed problem.

# Testbed problem

We take a standard SIR model:

$$\begin{aligned}
S' &= -\beta SI/N \\
I' &= \beta SI/N - \alpha I\\
R' &= \alpha I\\
N &= S + I + R\\
\end{aligned}$$

and solve it forward it time over $t \in [0, 20]$ for the parameters $\theta_{true}$ and initial conditions ${y_0}_{true}$ to get a true solution $y_{true}$:

$$\begin{aligned}
\theta_{true} = [\beta_{true}, \alpha_{true}] &= [1.3, 0.2]\\
{y_0}_{true} = [{S_0}_{true}, {I_0}_{true}, {R_0}_{true}] &= [999, 1, 0]
\end{aligned}$$

We generate (noisy) data $y$ by observing only the states $[S, R]$ over a time period $t \in [0, 7]$, which we encapsulate in an observation function, $g(\cdot)$, subject to Poisson noise:

$$y \sim \text{Poisson}(g(y_{true}))$$

# Standard variant

We fit the objective function $H_{standard} = -2\log\mathcal{L}_{standard}$:

$$
H_{standard} = \lVert L(y - g(\Phi c))\rVert^2 + \lVert W(D\Phi c - f(\Phi c, \theta)) \rVert^2 - 2\log{|L|} - 2\log{|W|}
$$
where
$$L = \frac{1}{\sigma_L} \mathbb{I},\qquad W = \frac{1}{\sigma_W} \mathbb{I}$$
and we define the weights:
$$w = \left[\frac{1}{\sigma_L}, \frac{1}{\sigma_W}\right]$$

We solve this iteratively:

1. Set the iteration number $i=0$.
2. Set initial guess $c^{(0)}, \theta^{(0)}$ and initial weights $w^{(0)}$.
3. Minimise $H_{standard}$ over $c, \theta$ for given weights $w^{(i)}$ $\to c^{(i+1)}, \theta^{(i+1)}$.
4. Minimise $H_{standard}$ over $w$ for given $c^{(i+1)}, \theta^{(i+1)}$ $\to w^{(i+1)}$.
5. If $i$ \< maximum iterations, set iteration number $i \gets i+1$ and return to step 2.

We shortcut step 4 by instead doing:

$w^{(i+1)} \gets s(c^{(i+1)}, \theta^{(i+1)})$
  
where $s(c, \theta)$ = $\sqrt{\frac{n}{r(c, \theta)}}$ and $r$ is the residual function, $n$ the length of $r(c, \theta)$.

In this case $r(c, \theta) = \left[y - g(\Phi c), \;D\Phi c - f(\Phi c, \theta)\right]$. We justify this by stating that this is the `optimal' solution via solution of the KKT conditions, once $r$ is known (assuming there are no constraints imposed on the problem, which is not necessarily true). 

Here, we `absorb' scaling into the covariance matrices, but because we do not optimise for them directly, we don't get the proper scalings.

## Results
We use initial weights $[1, 1]$, and initial guesses $\sim \text{Poisson}(1000)$.
We have to terminate the iteration early, at iteration $i=1$ where $w=[0.05, 0.49]$, and recover parameter estimates $\beta = 1.19, \alpha = 0.19$. 

![Fit to $H_{standard}$ at iteration $i=1$ ($\beta=1.19, \alpha=0.19$)](img/standard_fit.png)

The model weights seem to grow at each iteration. This eventually leads to underfits.

![Weights for each iteration, solving $H_{standard}$](img/standard_weights.png)

Iteration 0 seems slightly overfit, due to the oscillations in the estimate at small time in $S$. By iteration 2, the fit is already degraded, and by iteration 3, the estimate is very underfit.

![Fit to $H_{standard}$ at iteration $i=0$ ($\beta=1.21, \alpha=0.20$)](img/standard_fit_iter0.png)

![Fit to $H_{standard}$ at iteration $i=2$ ($\beta=1.00, \alpha=0.16$)](img/standard_fit_iter2.png)

![Fit to $H_{standard}$ at iteration $i=3$ ($\beta=0.70, \alpha=0.10$)](img/standard_fit_iter3.png)



# Delta-t variant
We fit the objective function $H_{\Delta t} = -2\log\mathcal{L}_{\Delta t}$:

$$
H_{standard} = \lVert L(y - g(\Phi c))\rVert^2 + \lVert W(\,\Delta t(D\Phi c - f(\Phi c, \theta))\,) \rVert^2 - 2\log{|L|} - 2\log{|W|}
$$
where
$$L = \frac{1}{\sigma_L} \mathbb{I},\qquad W = \frac{1}{\sigma_W} \mathbb{I}$$
and we define the weights:
$$w = \left[\frac{1}{\sigma_L}, \frac{1}{\sigma_W}\right]$$

Essentially, we introduce a scaling $\Delta t$ into the _residual_ function.
This fixes the consistency problem with the SDE interpretation.


