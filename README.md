# TFG_Physics
# GBM occupation time simulations

This repository contains the code used to simulate the occupation time of a geometric Brownian motion inside an interval.

The process is

```math
dX_t = \mu X_t dt + \sigma X_t dW_t,
```

and the quantity studied is the occupation time

```math
T_{[a,b]}(t) = \int_0^t \mathbf{1}_{[a,b]}(X_s)\,ds.
```

In words, this measures how much time the trajectory spends inside the interval `[a,b]`.

## Goal

The goal of the code is to compare Monte Carlo simulations with the theoretical long-time formulas for:

* the first moment,

```math
\mathbb{E}[T_{[a,b]}(t)],
```

* the second moment,

```math
\mathbb{E}[T_{[a,b]}(t)^2],
```

* and the ergodicity breaking parameter,

```math
EB = \frac{\mathbb{E}[T_{[a,b]}(t)^2]}{\mathbb{E}[T_{[a,b]}(t)]^2} - 1.
```

## Files

The main files are:

```text
gbm_occupation_simulation.ipynb
```

Notebook with the simulations, theoretical expressions and figures.

```text
figures/
```

Folder where the generated figures are saved.
## Simulation method

The geometric Brownian motion is simulated using its logarithmic form. The code simulates the process

```math
Y_t = \log X_t.
```

For a GBM,

```math
dX_t = \mu X_t dt + \sigma X_t dW_t,
```

the logarithm satisfies

```math
Y_{t+h}
=
Y_t
+
\left(\mu-\frac{\sigma^2}{2}\right)h
+
\sigma \sqrt{h} Z,
```

where (Z) is a standard normal random variable.

In the code, each trajectory starts from

```math
Y_0 = \log(x_0).
```

Then, at every time step, a new normal random variable is generated and the logarithmic process is updated. Since the occupation time is defined as the time spent by (X_t) inside the interval ([a,b]), the code checks the equivalent condition

```math
\log(a) \leq Y_t \leq \log(b).
```

If this condition is satisfied, the occupation time is increased by the time step (h). Therefore, the occupation time is approximated numerically as

```math
T_{[a,b]}(t)
\approx
\sum_k \mathbf{1}_{[\log(a),\log(b)]}(Y_{t_k})h.
```

This procedure is repeated for many independent Monte Carlo trajectories. For each value of (x_0) and (\mu), the code stores the simulated occupation times at several observation times. Then the first moment, second moment and ergodicity breaking parameter are estimated using the sample averages:

```math
\mathbb{E}[T_{[a,b]}(t)]
\approx
\frac{1}{N}\sum_{i=1}^N T_i(t),
```

```math
\mathbb{E}[T_{[a,b]}(t)^2]
\approx
\frac{1}{N}\sum_{i=1}^N T_i(t)^2,
```

and

```math
EB(t)
=
\frac{\mathbb{E}[T_{[a,b]}(t)^2]}
{\mathbb{E}[T_{[a,b]}(t)]^2}
-1.
```

## Parameters used

The main parameters are:

```python
a = 1.0
b = 2.0
sigma = 1.0
```

Different values of `mu` are used to study the three regimes:

```python
mu_values = [0.0, 0.5, 1.0]
```

which correspond to:

```math
\sigma^2 > 2\mu,
```

```math
\sigma^2 = 2\mu,
```

and

```math
\sigma^2 < 2\mu.
```

Different initial conditions are also used:

```python
x0_values = [0.7, 1.2, 2.3]
```

These correspond to the three possible regions:

* `x0 < a`,
* `a <= x0 <= b`,
* `x0 > b`.

## How to run

First install the required Python packages:

```bash
pip install -r requirements.txt
```

Then open the notebook:

```bash
jupyter notebook gbm_occupation_simulation.ipynb
```

or run it in Jupyter Lab:

```bash
jupyter lab
```

## Requirements

The code uses:

```text
numpy
matplotlib
numba
jupyter
```

A minimal `requirements.txt` is:

```text
numpy
matplotlib
numba
jupyter
```

## Output

The code produces:

* Monte Carlo estimates of the first moment,
* Monte Carlo estimates of the second moment,
* Monte Carlo estimates of the EB parameter,
* theoretical values,
* comparison plots between simulation and theory.

The figures are saved in the `figures/` folder.

## Notes

The time step `h` and the number of simulated paths `n_paths` control the accuracy of the simulation.

Smaller `h` gives a better approximation of the continuous-time process, but makes the simulation slower.

Larger `n_paths` gives smoother Monte Carlo averages, but also increases the computation time.
## Conditional PDF simulation

The code also computes the limiting probability density of the occupation time. In the transient regimes, the limiting distribution can have a discrete mass at zero. This mass corresponds to the probability that the GBM never reaches the interval ([a,b]).

Therefore, the limiting law can be written as

```math
P(T^\infty_{[a,b]} \in du)
=
P(\tau_I=\infty)\delta_0(du)
+
P(\tau_I<\infty)f_\infty(u\mid \tau_I<\infty)du.
```

The Monte Carlo simulation used for the density computes the conditional distribution (f_\infty(u\mid \tau_I<\infty)). For this reason, when the initial condition starts outside the interval, the simulated path is started at the entrance boundary: at (\log(a)) if (x_0<a), and at (\log(b)) if (x_0>b). This uses the Markov property and avoids simulating paths that never enter the interval.

The process is simulated in logarithmic variables,

```math
Y_t=\log X_t,
```

using

```math
Y_{t+\Delta t}
=
Y_t+
\left(\mu-\frac{\sigma^2}{2}\right)\Delta t
+
\sigma\sqrt{\Delta t}Z.
```

At each step, the code checks whether

```math
\log(a)\leq Y_t\leq \log(b).
```

If this condition is satisfied, (\Delta t) is added to the occupation time. The simulation is run for a large final time, and in the transient regimes the path is stopped once it has moved far enough in the direction of the effective drift, so that the probability of returning to the interval is negligible.

The theoretical density is computed as a series of exponentials,

```math
f_\infty(u\mid \tau_I<\infty)
=
\sum_n A_n e^{-\lambda_n u},
```

where the values (\lambda_n) are obtained from the spectral equation and the coefficients (A_n) are computed from the residues of the conditional Laplace transform.

