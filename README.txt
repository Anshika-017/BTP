Gamma-Ray Detector Calibration and Parameter Estimation

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

$$
x=R\cos(\theta), \qquad y=R\sin(\theta)
$$

The angular value is converted from degrees to radians during the calculation.

### 2. Calibration Point Weights

The experimentally measured counts \(M_i\) for the 108 calibration points are read from `det1.txt`. The relative weight of each point is:

$$
W_i=\frac{\sqrt{M_i}}{\sqrt{\sum_{j=1}^{108}M_j}}
$$

where:

- \(W_i\) is the weight assigned to the \(i\)-th calibration point.
- \(M_i\) is the experimentally measured count at the \(i\)-th point.

The weight represents the relative importance of each calibration point when the overall error is calculated.

### 3. Detector Reconstruction

The function `reconF()` in `reconF.h` calculates the theoretical detector count for a given calibration position and a selected set of fitting parameters.

It uses the detector geometry to determine the relevant distances and geometrical quantities before calling the Monte Carlo simulation.

### 4. Distance Travelled and Interaction Probability

For each simulated gamma ray, the program calculates the distance travelled inside the detector crystal. This path length is used to calculate the probability of interaction inside the detector:

$$
P=1-e^{-\mu d}
$$

where \(d\) is the calculated path length and \(\mu\) is the total attenuation coefficient.

### 5. Monte Carlo Simulation and Efficiency

The Monte Carlo model estimates the detector efficiency by simulating possible gamma-ray trajectories and their interactions inside the detector.

The simulation considers:

- Detector geometry and path lengths.
- Angular sampling.
- Gamma-ray attenuation.
- Photoelectric absorption and Compton interaction contributions.
- Attenuation through the detector wall.

The simulation uses 12,000 Monte Carlo histories per calculation in the current implementation.

For each simulated gamma ray, the program:

1. Calculates the interaction probability inside the detector.
2. Accounts for attenuation through the detector wall.
3. Calculates the angular weights (theta weight \(w_\theta\) and alpha weight \(w_\alpha\)).
4. Multiplies the interaction probability by both weights and accumulates the result.

The final efficiency is:

$$
\varepsilon=\frac{1}{N}\sum_{i=1}^{N}P_i\,w_{\theta,i}\,w_{\alpha,i}
$$

where:

- \(N=12{,}000\) is the number of Monte Carlo iterations.
- \(P_i\) is the interaction probability of the \(i\)-th simulated gamma ray.
- \(w_{\theta,i}\) is the theta weight (angular weight for the polar angle) of the \(i\)-th gamma ray.
- \(w_{\alpha,i}\) is the alpha weight (angular weight for the angle \(\alpha\)) of the \(i\)-th gamma ray.

\(P_i\) includes the detector interaction and the attenuation through the detector wall:

$$
P_i=\left(1-e^{-\mu d_i}\right)e^{-\mu_{\text{wall}}t_i/\cos\theta_i}
$$

where \(d_i\) is the path length travelled inside the detector, \(\mu\) is the total attenuation coefficient of the detector, \(\mu_{\text{wall}}\) is the attenuation coefficient of the wall, \(t_i\) is the wall thickness, and \(\theta_i\) is the angle of the ray.

The final efficiency is stored in the variable `inteff` in the program, and is used to calculate the theoretical detector count.

### 6. Parameter Estimation

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

### 7. Error Calculation

**Relative error term.** For each of the 108 calibration points, the symmetric relative difference between the theoretical and experimental counts is calculated:

$$
t_i=\frac{C_i-M_i}{(C_i+M_i)/2}
$$

where:

- \(C_i\) is the theoretically predicted count.
- \(M_i\) is the experimentally measured count.

**Weighted error contribution.** The relative difference is squared and multiplied by the corresponding weight:

$$
\phi_i=W_i\,(t_i)^2
$$

**Total error.** The individual contributions are summed over all 108 calibration points:

$$
E_{\text{total}}=\sum_{i=1}^{108}\phi_i
$$

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

$$
E_{\text{best}}\approx 1.61
$$

The corresponding parameter values are approximately:

- `ss = 0.00055`
- `tau = 3 × 10⁻⁶`
- `wallmu = 0.0986856`

The desired error threshold is:

$$
E_{\text{target}}=0.006
$$

The current implementation has not yet reached this target. The reported result is the best value found within the tested parameter search and trial limit, not a converged solution.

## Future Work

### 1. Optimisation of the Error Function

- Investigate alternative error calculation methods and different combinations of error terms.
- Compare their performance to identify an error function that minimises the difference between theoretically predicted and experimentally measured counts.
- Evaluate the effect of different weighting approaches on the overall error.

### 2. Improvement of Random Number Selection and Sampling

- Experiment with alternative random-number selection strategies for choosing the fitting parameters `wallmu`, `ss`, and `tau`.
- Explore different sampling techniques instead of relying only on the current random selection within fixed parameter ranges.
- Investigate modifications to random-number generation and sampling in the Monte Carlo simulation to better represent possible gamma-ray trajectories.
- Assess whether improved parameter sampling and Monte Carlo sampling lead to more stable efficiency estimates and lower total error.

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
