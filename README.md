# GARCH-J-model
This project implements and analyzes a GARCH-J model using both simulated and real data (log returns of the OMX index).

Model Overview
The GARCH-J model extends the standard GARCH model by including jumps in returns.

Steps
Simulate Data: Generate 1000 data points using specified parameters for jump probability, mean return, and volatility.

Estimate Parameters:

QMLE: Estimate model parameters using a Gaussian approximation.
Exact MLE: Estimate model parameters using exact maximum likelihood.
Monte Carlo Simulation: Perform 100 repetitions to analyze the distribution of estimated parameters.

Fit Real Data: Apply both QMLE and exact MLE methods to the real OMX log return data.
