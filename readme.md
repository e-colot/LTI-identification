# Goal

The goal of this lab is for you to gain experience with the tools for identifying LTI systems described during the lectures. This will be applied to measurements on a real dynamic system. Either, you can use one of the systems available at the department, or you can suggest another system, e.g. which would be of interest for your master thesis.

# Typical project

The tasks enumerated below are a minimal set of tasks and points of interest to be considered in this lab:

1. Design an excitation signal, for instance a multisine with holes, to trap nonlinear distortions,

2. Apply the signal to the system. Does the response match your expectations? Possibly, adapt the excitation signal (RMS, excited frequency band, frequency resolution, choice of non-excited frequencies, …)

3. Distinguish the nonlinear distortions from the noise via the use of the robust and/or the fast method,

4. Obtain a non-parametric estimate of the BLA of the system, with a quantification of its uncertainty

5. Implement an identification procedure to obtain a parametric model of the BLA.

6. Apply model order selection criteria (AIC, whiteness test of the residuals, value of the cost function,…). Consider the bias-variance trade-off.

7. Be critical of your own results. Does the estimated parametric model satisfy your expectations?

Additional points of interest are

- make experiments at different- setpoints, and compare the BLAs- (important if you want to- classify the nonlinear systems,- see the course on nonlinear- systems),
- perform an experiment where the- setpoint is varying during the- experiment, and compute the- BLTIA (Best Linear Time- Invariant Approximation) of the- system.
- Determine and validate the- uncertainty of the estimated- parametric model
- Open loop versus closed loop- experiments
- Band limited versus- zero-order-hold assumption
- ZOH: model from reference to- output, without anti-alias- filter
- BL: model from measured input to measured output, with AA filter

# Expectations

Prepare a written, well documented report of your measurement campaign and your identification results.
Possible systems to consider

## Wiener-Hammerstein

This is a system consisting of two linear blocks, with a static nonlinearity in between. (Keep the amplitude of the input signal smaller than 1.5V.)

## Silverbox

This is an electronic circuit which emulates a mass-spring damper system, with a nonlinear spring.

## Linear parameter-varying system

which has an additional input which allows to vary one of the physical parameters of the system.

