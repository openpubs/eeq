
# EEQ Model Simulation for Exposure Selection in Survival Analysis

This project performs a simulation study using an Empirical Exposure Quantile (EEQ) model in the context of survival analysis. The goal is to estimate the contribution of multiple exposures to a survival outcome and identify key exposures based on weighted indices.

## Overview

The script simulates survival data with covariates and applies an optimization-based method to estimate:
- A coefficient for an exposure index (`EEQ`)
- Weights for each individual exposure contributing to the index
- Whether each exposure is "key", based on a weight threshold

The optimization is based on maximizing a penalized Cox partial likelihood with a constraint that exposure weights sum to 1.

## Features

- Simulates multivariate exposure data with defined correlation
- Introduces additional covariates (`Z1`, `Z2`)
- Generates right-censored survival times
- Computes empirical quantiles of exposures
- Optimizes an EEQ-based Cox proportional hazards model
- Identifies and tracks key exposures across multiple simulations

## File Structure

- **`EEQ.R`**: Main script containing the entire simulation and estimation pipeline.
- **Outputs**: 
  - `result.output`: Proportion of simulations where each exposure was identified as "key".

## Main Components

### Data Generation

- 10 continuous exposures (`X1`–`X10`) from a multivariate normal distribution
- Binary covariate `Z1` and continuous covariate `Z2`
- True coefficient vector assigns signal to the first 3 exposures and `Z1`, `Z2`

### Survival Time Simulation

- Survival times generated from exponential distribution with hazard dependent on exposures and covariates
- Censoring introduced via exponential distribution

### EEQ Model

- Constructs weighted exposure index from empirical quantiles of exposures
- Fits a penalized Cox proportional hazards model
- Penalizes deviation of weight sum from 1
- Uses `nlminb()` for optimization

### Output

- Estimated EEQ coefficients across simulations
- Exposure weights per simulation
- Frequency each exposure is selected as “key” (weight > 1/p)

## Requirements

- Packages:
  - `survival`
  - `MASS`
  - `glmnet`

Install missing packages using:

```r
install.packages(c("survival", "MASS", "glmnet"))
```

## How to Run

1. Open `EEQ.R` in RStudio or any R environment.
2. Set your desired simulation parameters:
   - `n`: number of observations (default = 1000)
   - `p`: number of exposures (default = 10)
   - `Nsim`: number of simulation iterations (default = 3)
3. Run the script.
4. The final output, `result.output`, shows the proportion of times each exposure was selected.

## Example Output

```
      eeq.select.weight
[1,]              1.000
[2,]              1.000
[3,]              1.000
[4,]              0.000
[5,]              0.000
...
```

## License

This project is intended for academic and research purposes. Please cite appropriately if used in publications.
