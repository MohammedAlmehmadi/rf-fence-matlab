# RF Fence — MATLAB Physics Simulation

**A radio-frequency fence simulation framework for electromagnetic diffraction modeling and signal-detection research.**

## Overview

RF Fence investigates how an opaque object crossing the line of sight between a transmitter and receiver affects the received electromagnetic field.

This repository contains the physics layer of the project, including diffraction models, numerical integration, coordinate handling, and a staged validation framework.

**Project status:** Physics verification in progress. Expected numerical values were obtained through a NumPy transcription of the MATLAB implementation. Execution in MATLAB and experimental validation remain pending.

## Key Features

- Scalar diffraction modeling using the Babinet/Fresnel paraxial approximation.
- Numerical Kirchhoff integration for comparison with approximate models.
- Toolbox-free Fresnel C/S integral calculations.
- A unified scattered-field interface using line-of-sight-relative coordinates.
- An axial opaque-disk reference model for diffraction analysis.
- Staged numerical checks covering units, convergence, symmetry, and model consistency.

## Repository Structure

| File | Description |
|---|---|
| `fsr_wavelength.m` | Calculates wavelength from operating frequency. |
| `fsr_fresnel_cs.m` | Computes Fresnel C/S integrals. |
| `fsr_babinet_paraxial.m` | Implements the paraxial Babinet/Fresnel model. |
| `fsr_kirchhoff_numeric.m` | Implements numerical Kirchhoff integration. |
| `fsr_scattered_field.m` | Provides a unified scattered-field calculation interface. |
| `fsr_disk_axial.m` | Models the axial diffraction field of an opaque disk. |
| `fsr_validate_physics.m` | Runs staged physics validation checks. |

## Getting Started

### Requirements

- MATLAB
- All seven `.m` files in the same directory

### Run the validation script

1. Clone or download this repository.
2. Open MATLAB and set the Current Folder to the downloaded directory.
3. Run:

```matlab
fsr_validate_physics
```

4. Review the validation output and compare it with the reference values in `README_MATLAB_v2.txt`.

## Validation Approach

The validation framework is organized into stages covering:

- Units and geometric parameters.
- Numerical convergence and limiting cases.
- Symmetry, reciprocity, and coordinate consistency.
- Paraxial versus numerical diffraction models.
- Obstruction-size sensitivity.
- Approximation accuracy across different observation positions.
- Sensitivity to operating frequency and geometry.

Reference values are documented in the accompanying text file. These values have not yet been confirmed by executing the MATLAB implementation.

## Scope and Limitations

This is a scalar physical-optics diffraction model for an opaque obstruction. It is not a full-wave electromagnetic solver.

The current implementation has not been validated against laboratory measurements or field experiments. Real-world detection performance, false-alarm rates, and operational reliability remain to be established.

## Roadmap

- [ ] Execute and verify the physics validation suite in MATLAB.
- [ ] Resolve discrepancies and document verified test results.
- [ ] Validate and compare signal-detection algorithms.
- [ ] Generate reproducible plots and numerical performance metrics.
- [ ] Integrate software-defined radio (SDR) hardware.
- [ ] Evaluate the model against real RF measurements.

## Project Status

**Current phase: Physics-layer verification.**

The goal is to establish a validated simulation foundation before proceeding to hardware integration and experimental testing.
