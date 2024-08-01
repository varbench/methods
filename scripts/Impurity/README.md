# Quantum Impurity Model Scripts

This repository contains scripts for calculating the ground state energy and variance of various quantum impurity models using forkTPS. 

## Overview

The scripts in this repository perform the following steps:
1. Load the corresponding hybridization function for the chosen model.
2. Perform a deterministic discretization to obtain the impurity Hamiltonian parameters.
3. Calculate the ground state energy and variance using DMRG.

## File Naming Convention

The scripts follow a consistent naming convention: `model_Nb.md`

Where:
- `model` represents the name of the quantum impurity model
- `Nb` indicates the number of bath sites per spin-orbital

## Available Models and Hybridization Functions

1. **Single-Band Model**
   - Hybridization function stored in: `OneBand_InFile`

2. **Three-Band Model without Spin-Orbital Coupling**
   - Hybridization function stored in: `SrRuO214_NoSOC_InFile`

3. **Three-Band Model with Spin-Orbital Coupling**
   - Hybridization function stored in: `SrRuO214_SOC_InFile`