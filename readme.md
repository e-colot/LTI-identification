# Goal

The goal of this lab is for you to gain experience with the tools for identifying LTI systems described during the lectures. This will be applied to measurements on a real dynamic system. Either, you can use one of the systems available at the department, or you can suggest another system, e.g. which would be of interest for your master thesis.

Prepare a written, well documented report of your measurement campaign and your identification results.

# Typical project

> The tasks enumerated below are a minimal set of tasks and points of interest to be considered in this lab:
> 
> 1. Design an excitation signal, for instance a multisine with holes, to trap nonlinear distortions,
> 
> 2. Apply the signal to the system. Does the response match your expectations? Possibly, adapt the excitation signal (RMS, excited frequency band, frequency resolution, choice of non-excited frequencies, …)
> 
> 3. Distinguish the nonlinear distortions from the noise via the use of the robust and/or the fast method,
> 
> 4. Obtain a non-parametric estimate of the BLA of the system, with a quantification of its uncertainty
> 
> 5. Implement an identification procedure to obtain a parametric model of the BLA.
> 
> 6. Apply model order selection criteria (AIC, whiteness test of the residuals, value of the cost function,…). Consider the bias-variance trade-off.
> 
> 7. Be critical of your own results. Does the estimated parametric model satisfy your expectations?

> Additional points of interest are
> 
> - make experiments at different- setpoints, and compare the BLAs- (important if you want to- classify the nonlinear systems,- see the course on nonlinear- systems),
> - perform an experiment where the- setpoint is varying during the- experiment, and compute the- BLTIA (Best Linear Time- Invariant Approximation) of the- system.
> - Determine and validate the- uncertainty of the estimated- parametric model
> - Open loop versus closed loop- experiments
> - Band limited versus- zero-order-hold assumption
> - ZOH: model from reference to- output, without anti-alias- filter
> - BL: model from measured input to measured output, with AA filter


# How to use the PXI measurement system

1. If necessary, turn on the PXI system (the computer turns on automatically), and the supply voltage of the system under test (if applicable).
2. Open the VI in Labview
    > 1. Open the Program YR_Run_IF.vi labview VI (Virtual Instrument)
    > 2. Run the VI
3. Set the signal directory to your own dir.
4. Prepare the signals
    * Either use the built in signal generator (multisine)
    * Or prepare your own signal
        > The excitation signal needs to be prepared as a matrix, where each column represents a realisation. Thus, the matrix has dimension Npp x R with Npp the number of samples (points) per period, and R the number of realisations.
        >
        >> - Technical note: only odd columns will be considered as realisations. So, you should put a ‘dummy’ column at each even index.
        >> - Technical note: at least 2 realisations must be provided for the system to work properly
        >
        > Save the input signal matrix in a file named ***[filename]_Sig_E0_S0.mat*** where [filename] can be chosen freely.
        > You also need to prepare a matrix with the excited frequencies (expressed in bins)
        > Save the vector of excited bins in a file named ***[filename]_Sel_E0_S0.mat***

    * Set the sampling frequency and the number of periods in the VI.
    * Fill in [filename] in the field “Generic Signal Filename”
5. Load the signal: Button Calc Excitations
6. Set the power sweep, and load it: Button Load
7. Initialise the device: Button Initialize
8. Measure
    > Plots of measurements are available
    > Signals are saved as Matlab files with this naming convention
    >>    ACQ_Rx_Py_Ez_M0_F0.mat: measured signal of Realisation x, Power y, E?, M?, F?
    >>    AWG_Rx_Ey.mat: Reference signal in AWG
    >>    AWGEXC_Rx_Ry.mat: Excited frequencies and detection lines (expressed in bin!)
    >>    REFACQ_Rx_Py_Ez_M0_F0.mat: extra information what is this?
    >
    >if necessary, adapt the range of the ACQ cards (tab ACQ) to the maximum value of the signal
9. To perform measurements with different settings, the VI needs to be restarted: stop and start the VI in LabView.


