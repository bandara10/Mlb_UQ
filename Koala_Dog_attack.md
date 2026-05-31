
# Statistical Analysis

## Observation Model

For each koala mortality record i:

Yi = 1 if death was due to the target cause

Yi = 0 otherwise

Assume:

Yi ~ Bernoulli(pi_i)

where:

pi_i = P(Yi = 1)

--------------------------------------------------

## Linear Predictor

eta_i = beta_0
      + beta_1*x_1i
      + beta_2*x_2i
      + ...
      + beta_p*x_pi
      + W(s_i)

where:

beta_0 = intercept

beta_j = regression coefficients

x_ij = explanatory variables

W(s_i) = spatial random effect

--------------------------------------------------

## Link Function

logit(pi_i)

= log[ pi_i / (1 - pi_i) ]

= eta_i

Equivalent form:

pi_i = exp(eta_i) / (1 + exp(eta_i))

--------------------------------------------------

## Bernoulli Likelihood

L(beta,W)

= Product over i

[pi_i^(y_i)] * [(1-pi_i)^(1-y_i)]

--------------------------------------------------

## Spatial Process

W(s)

~ Gaussian Process(0, C(theta))

where:

C(theta) = spatial covariance function

theta = spatial parameters

--------------------------------------------------

## Joint Likelihood

p(Y,W | beta, theta)

= p(Y | W, beta)

× p(W | theta)

--------------------------------------------------

## Prior Distributions

beta ~ p(beta)

theta ~ p(theta)

--------------------------------------------------

## Posterior Distribution

p(beta,W,theta | Y)

proportional to

p(Y | W,beta)

× p(W | theta)

× p(beta)

× p(theta)

--------------------------------------------------

## Computational Inference

Posterior inference performed using:

- INLA
- SPDE

--------------------------------------------------

MODEL HIERARCHY

Observed Data (Y)
        |
        v
 Bernoulli Model
        |
        v
    Likelihood
        |
        v
 Logistic Model
        |
        v
 Spatial Random Effect
        |
        v
 Joint Likelihood
        |
        v
      Priors
        |
        v
    Posterior
        |
        v
    INLA-SPDE
        |
        v
 Parameter Estimates



# Statistical Analysis

## 1. Observation Model

For each koala mortality record \(i\), define:

\[
Y_i =
\begin{cases}
1, & \text{if death was due to the target cause} \\
0, & \text{otherwise}
\end{cases}
\]

Assume:

\[
Y_i \sim \text{Bernoulli}(\pi_i)
\]

where:

\[
\pi_i = P(Y_i = 1)
\]

is the probability that mortality was due to the target cause.

---

## 2. Linear Predictor

The probability of cause-specific mortality was modelled using a logistic regression:

\[
\eta_i
=
\beta_0
+
\sum_{j=1}^{p}\beta_j x_{ij}
+
W(s_i)
\]

where:

- \(\eta_i\) = linear predictor
- \(\beta_0\) = intercept
- \(\beta_j\) = regression coefficients
- \(x_{ij}\) = explanatory variables
- \(W(s_i)\) = spatial random effect at location \(s_i\)

---

## 3. Link Function

The linear predictor was linked to the mortality probability using the logit link:

\[
\text{logit}(\pi_i)
=
\log
\left(
\frac{\pi_i}
{1-\pi_i}
\right)
=
\eta_i
\]

or equivalently:

\[
\pi_i
=
\frac{\exp(\eta_i)}
{1+\exp(\eta_i)}
\]

---

## 4. Bernoulli Likelihood

Conditional on the spatial effect:

\[
L(\boldsymbol{\beta},W)
=
\prod_{i=1}^{n}
\pi_i^{y_i}
(1-\pi_i)^{1-y_i}
\]

This likelihood measures how well the regression coefficients and spatial effects explain the observed mortality outcomes.

---

## 5. Spatial Process Model

Spatial dependence was modelled using a Gaussian random field:

\[
W(s)
\sim
GP(0,C(\theta))
\]

where:

- \(GP\) = Gaussian Process
- \(C(\theta)\) = covariance function
- \(\theta\) = spatial hyperparameters

---

## 6. Joint Likelihood

Combining the observation model and spatial process:

\[
p(\mathbf{Y},W|\boldsymbol{\beta},\theta)
=
p(\mathbf{Y}|W,\boldsymbol{\beta})
\times
p(W|\theta)
\]

This represents the joint likelihood of the observed mortality data and the latent spatial field.

---

## 7. Prior Distributions

Bayesian prior distributions were assigned:

\[
\boldsymbol{\beta}
\sim
p(\boldsymbol{\beta})
\]

\[
\theta
\sim
p(\theta)
\]

where:

- \(\boldsymbol{\beta}\) = regression parameters
- \(\theta\) = spatial covariance parameters

---

## 8. Posterior Distribution

The posterior distribution was obtained using Bayes' theorem:

\[
p(\boldsymbol{\beta},W,\theta|\mathbf{Y})
\propto
p(\mathbf{Y}|W,\boldsymbol{\beta})
\times
p(W|\theta)
\times
p(\boldsymbol{\beta})
\times
p(\theta)
\]

---

## 9. Computational Inference

Posterior inference was performed using:

- Integrated Nested Laplace Approximation (INLA)
- Stochastic Partial Differential Equation (SPDE) representation of the Gaussian random field

which provides an efficient approximation to the posterior distribution without requiring Markov chain Monte Carlo (MCMC) sampling.

---

# Model Hierarchy

```text
Observed mortality outcomes (Y)
                ↓
      Bernoulli distribution
                ↓
          Likelihood
                ↓
      Logistic regression
                ↓
     Spatial Gaussian field
                ↓
         Joint likelihood
                ↓
       Prior distributions
                ↓
     Posterior distribution
                ↓
          INLA-SPDE
                ↓
       Parameter estimates
```