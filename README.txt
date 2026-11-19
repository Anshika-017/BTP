# Gamma-Ray Detector Calibration and Parameter Estimation

## Overview

This project focuses on the computational modelling, reconstruction, and parameter estimation of a gamma-ray detector system.

The implementation combines a detector-response model with Monte Carlo simulation to estimate the theoretical detector counts at different calibration positions. These predicted counts are compared with experimentally measured counts, and an error function is used to evaluate the agreement between the two.

The program performs parameter estimation by randomly sampling selected model parameters and identifying the combination that produces the minimum total error.

## Objectives

- Model the response of a gamma-ray detector using its geometry and interaction probabilities.
- Calculate theoretical detector counts at specified calibration positions.
- Implement Monte Carlo simulation to estimate detector efficiency.
- Compare theoretical counts with experimentally measured detector counts.
- Estimate fitting parameters by minimising the total error.
- Study the influence of parameter selection and random sampling on the accuracy of the reconstruction.

## Repository Structure

```text
BTP/
│
├── calcu.cpp
│   └── Main parameter-estimation program.
│
├── reconF.h
│   └── Detector reconstruction and Monte Carlo
│       simulation functions.
│
├── crysd1.txt
│   └── Detector geometry and coordinate information.
│
├── det1.txt
│   └── Experimentally measured detector counts
│       at 108 calibration points.
│
├── gridp.txt
│   └── Calibration-point positions in cylindrical
│       coordinates (R, theta, z).
│
├── ss1_det1.txt
│   └── Sampled values of the parameter ss.
│
├── tau1_det1.txt
│   └── Sampled values of the parameter tau.
│
├── wallmu1_det1.txt
│   └── Sampled values of the parameter wallmu.
│
├── error1_det1.txt
│   └── Total error calculated for each trial.
│
└── README.md
    └── Project documentation.
```

## Methodology

### 1. Detector Geometry and Calibration Positions

The detector geometry is specified in `crysd1.txt`.

The calibration positions are provided in `gridp.txt` in cylindrical coordinates:

- R: Radial distance
- theta: Angular position in degrees
- z: Axial coordinate

These positions are converted into Cartesian coordinates before reconstruction:

\[
x=R\cos(\theta), \qquad y=R\sin(\theta)
\]

The angular value is converted from degrees to radians during the calculation.

### 2. Detector Reconstruction

The function `reconF()` in `reconF.h` calculates the theoretical detector count for a given calibration position and a selected set of fitting parameters.

It uses the detector geometry to determine the relevant distances and geometrical quantities before calling the Monte Carlo simulation.

### 3. Monte Carlo Simulation

The Monte Carlo model estimates the detector efficiency by simulating possible gamma-ray trajectories and their interactions inside the detector.

The simulation considers:

- Detector geometry and path lengths.
- Angular sampling.
- Gamma-ray attenuation.
- Photoelectric absorption and Compton interaction contributions.
- Attenuation through the detector wall.

The simulation uses 12,000 Monte Carlo histories per calculation in the current implementation.

The resulting efficiency is used to calculate the theoretical detector count.

### 4. Parameter Estimation

The main program, `calcu.cpp`, estimates the fitting parameters:

- `ss`
- `tau`
- `wallmu`

The current parameter ranges are:

| Parameter | Minimum | Maximum |
|---|---:|---:|
| wallmu | 0.0986856 | 0.0986856 |
| tau | 1.0 × 10⁻⁶ | 6.0 × 10⁻⁶ |
| ss | 0.00050 | 0.00060 |

The program randomly selects parameter values within the specified ranges and evaluates the detector response for all 108 calibration positions.

**Note:** `wallmu` is fixed in the current implementation because its minimum and maximum values are identical.

### 5. Error Calculation

For each calibration position, the predicted detector count is compared with the experimentally measured count.

The relative error is calculated as:

\[
t_i=\frac{C_i-M_i}{(C_i+M_i)/2}
\]

where:

- \(C_i\) is the theoretically calculated count.
- \(M_i\) is the experimentally measured count.

The squared relative error is combined with the corresponding calibration-point weight to obtain the contribution from each point.

The total error is obtained by summing the contributions from all 108 calibration points.

The parameter combination producing the lowest total error is retained as the best result.

## Requirements

- C++ compiler (G++)
- Linux or Windows Subsystem for Linux (WSL)
- Standard C++ libraries

No external third-party libraries are required by the current source code.

## Compilation and Execution

### Step 1: Clone the repository

```bash
git clone https://github.com/Anshika-017/BTP.git
cd BTP
```

### Step 2: Compile the program

```bash
g++ -O2 -o calcu calcu.cpp
```

### Step 3: Run the program

```bash
./calcu
```

The program performs a maximum of 2,000 trials by default.

To specify a different maximum number of trials, use the `MAXQ` environment variable:

```bash
MAXQ=300 ./calcu
```

For example, the above command runs a maximum of 300 parameter-estimation trials.

**Important:** Run the executable from the repository directory so that the program can access the input files.

## Output Files

The program generates the following files:

| File | Description |
|---|---|
| `ss1_det1.txt` | Sampled values of ss |
| `tau1_det1.txt` | Sampled values of tau |
| `wallmu1_det1.txt` | Sampled values of wallmu |
| `error1_det1.txt` | Total error corresponding to each trial |

The terminal also displays the trial number, updates when a new best error is found, and prints the best parameter combination at the end of execution.

## Current Results

In the current implementation, the best error obtained is approximately:

\[
E_{\text{best}}\approx 1.61
\]

The corresponding parameter values are approximately:

- `ss = 0.00055`
- `tau = 3 × 10⁻⁶`
- `wallmu = 0.0986856`

The desired error threshold is:

\[
E_{\text{target}}=0.006
\]

The current implementation has not yet reached this target. The reported result is the best value found within the tested parameter search and trial limit, not a converged solution.

## Future Work

### 1. Optimisation of the Error Function

- Explore alternative combinations of error calculation terms.
- Investigate different approaches to comparing theoretical and experimental counts.
- Identify an error function that provides better agreement and reduces the total error.

### 2. Improvement of Random Number Selection

- Experiment with alternative random-number selection methods during parameter estimation.
- Explore improved sampling strategies for selecting fitting parameters.
- Investigate modifications to random-number generation in the Monte Carlo simulation to improve sampling and potentially reduce the total error.

## Notes

- The input files must be present in the working directory before execution.
- The program currently evaluates 108 calibration positions.
- The maximum number of parameter-estimation trials can be controlled using `MAXQ`.
- The detector geometry and parameter search ranges are specific to the current implementation.
- The numerical results depend on the input data, detector geometry, parameter ranges, and Monte Carlo sampling.

## Author

**Anshika Singh**

B.Tech – Chemical Engineering  
Indian Institute of Technology Jammu

**Project:** B.Tech Project (BTP)
