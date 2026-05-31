STATISTICAL MODEL STRUCTURE

Step 1. Define the outcome

For each koala:

Yi = 1 if death was due to dog attack

Yi = 0 if death was due to another cause

(or vehicle collision versus all other causes)

--------------------------------------------------

Step 2. Specify the probability model

Assume each observation follows a Bernoulli distribution.

This means each koala death can have only two outcomes:

- Case
- Control

Mathematically:

Yi ~ Bernoulli(pi)

where pi is the probability that death was due to the target cause.

--------------------------------------------------

Step 3. Relate probability to explanatory variables

The probability is not assumed to be constant.

Instead:

Probability depends on

- Dog density
- Road density
- Human population density
- Elevation
- Other risk factors

using:

logit(pi)

=
intercept
+
effects of explanatory variables
+
spatial random effect

--------------------------------------------------

Step 4. Construct the likelihood function

Given the Bernoulli probability model,

the likelihood for all observations is:

Likelihood

=
Probability of observing all recorded deaths

given the model parameters.

This is called the Bernoulli likelihood.

The likelihood tells us:

"How well do the current parameter values explain the observed mortality data?"

--------------------------------------------------

Step 5. Add spatial dependence

Nearby koala deaths may be more similar than distant deaths.

To account for this:

A spatial Gaussian random field W(s) is added.

This captures unexplained spatial clustering.

Examples:

- Local habitat characteristics
- Unmeasured road hazards
- Local dog activity

--------------------------------------------------

Step 6. Construct the joint likelihood

The model now contains two probabilistic components:

1. Bernoulli mortality process

2. Spatial Gaussian process

These are combined into a single joint likelihood.

Joint likelihood

=
Mortality likelihood

×

Spatial process likelihood

--------------------------------------------------

Step 7. Add prior distributions

Because the analysis is Bayesian,

prior distributions are assigned to:

- Regression coefficients
- Spatial parameters

These priors represent knowledge before seeing the data.

--------------------------------------------------

Step 8. Obtain the posterior distribution

Bayes' theorem combines:

Likelihood

×

Prior information

to produce:

Posterior distribution

The posterior contains everything known about the parameters after observing the data.

--------------------------------------------------

Step 9. Estimate parameters

The posterior distribution is approximated using:

INLA-SPDE

This avoids computationally intensive MCMC simulation.

--------------------------------------------------

COMPLETE MODEL JOURNEY

Observed koala deaths
          ↓
Bernoulli probability model
          ↓
Bernoulli likelihood
          ↓
Add explanatory variables
          ↓
Logistic regression model
          ↓
Add spatial Gaussian process
          ↓
Joint likelihood
          ↓
Add priors
          ↓
Posterior distribution
          ↓
INLA-SPDE computation
          ↓
Odds ratios, credible intervals,
spatial risk maps and predictions