# OPENFOAM TWO-TEMPERATURE THERMOCHEMICAL NONEQUILIBRIUM MODEL

A comprehensive OpenFOAM implementation of a two-temperature thermochemical non-equilibrium model for high-enthalpy air, designed for hypersonic reacting flow simulations.

## Table of Contents
- [Overview](#overview)
- [Governing Equations](#governing-equations)
- [Key Features](#key-features)
- [Project Structure](#project-structure)
- [Getting Started](#getting-started)
- [Implementation Details](#implementation-details)
- [Results](#results)
- [References](#references)

## Overview

This project implements a **two-temperature thermochemical non-equilibrium model** within the OpenFOAM framework, based on the seminal work by Casseau et al. (2016). The model is specifically designed for high-enthalpy air applications in hypersonic flows where vibrational and translational temperatures are not in equilibrium.

### Objectives
- Analyze one-temperature Navier–Stokes formulation baseline
- Develop two-temperature energy decomposition
- Implement physics within OpenFOAM using Mutation++
- Validate against literature results

## Governing Equations

### 1. Navier–Stokes Equations

**Mass conservation:**
```
∂ρ/∂t + ∇·(ρu) = 0
```

**Momentum conservation:**
```
∂(ρu)/∂t + ∇·(ρu⊗u) = −∇p + ∇·τ + ρf
```

**Energy conservation:**
```
∂(ρE)/∂t + ∇·[u(ρE + p)] = ∇·(k∇T) + ∇·(τ·u) + ρf·u
```

**Total energy:**
```
E = e + 1/2 |u|²
```

**Equation of state (ideal gas):**
```
p = ρRT
```

### 2. Multi-Species Formulation

Air is treated as a multi-species mixture (5 species: **air5 model**).

**Species transport equation:**
```
∂(ρYᵢ)/∂t + ∇·(ρuYᵢ) = −∇·Jᵢ + ω̇ᵢ
```

**Constraint:**
```
Σ Yᵢ = 1
```
Species mass fractions are strongly coupled with temperature and density:
```
Yᵢ ↔ T ↔ ρ ↔ p ↔ u
```

### 3. Two-Temperature Model

The fundamental innovation: **internal energy decomposition**

**Energy split:**
```
e(T, Tᵥ, Y) = eₜ(T,Y) + eᵥ(Tᵥ,Y)
```
Where:
- **eₜ** = translational–rotational energy (function of translational temperature T)
- **eᵥ** = vibrational–electronic energy (function of vibrational temperature Tᵥ)

**Total energy:**
```
E = eₜ + eᵥ + 1/2 |u|²
```

**System of equations:**
1. Conservative total energy equation
2. Vibrational energy equation:
```
deᵥ/dt = Qₜᵣ↔ᵥ + Qcₕₑₘ→ᵥ
```
This accounts for:
- **Vibration–translation relaxation** (energy exchange between modes)
- **Chemical energy exchange** (coupling with reactions)

## Key Features

✅ **Two-temperature energy decomposition** – Separate translational and vibrational temperatures  
✅ **Multi-species support** – 5-species air model with chemical reactions  
✅ **Mutation++ integration** – Advanced thermochemical nonequilibrium library  
✅ **OpenFOAM compatibility** – Built on ShockThermo solver  
✅ **Parallel computing** – MPI support with effective speedup up to 4 processes  
✅ **Validated results** – Excellent agreement with reference literature  

## Project Structure

```
.
├── applications/      # Solver applications and utilities
├── bin/              # Compiled executables
├── etc/              # Configuration and settings
├── src/              # Source code
│   ├── thermophysicalModels/    # Two-temperature thermophysics
│   ├── solvers/                 # Custom solver implementations
│   └── ...
├── thirdParty/       # External libraries (Mutation++)
├── tutorials/        # Example cases and test scenarios
├── test/            # Test cases and validation
├── Allwmake         # Build script
└── README.md        # This file
```

## Getting Started

### Prerequisites
- OpenFOAM (compatible version)
- Mutation++ library
- C++ compiler (g++, clang, or ICC)
- MPI library (OpenMPI or MPICH)
- Python 3 (optional, for post-processing)

### Building

```bash
# Source OpenFOAM environment
source $FOAM_INSTALL_DIR/etc/bashrc

# Build the project
./Allwmake
```

### Running a Test Case

```bash
# Navigate to tutorials
cd tutorials/singleCell

# Prepare case
blockMesh
# or use existing mesh

# Run solver (single process)
shockThermoTwo

# Or run in parallel
mpirun -np 4 shockThermoTwo -parallel
```

## Implementation Details

### OpenFOAM Modifications

**Base solver:** ShockThermo (derived from ShockFluid)

**Key function:** `thermophysicalPredictor()`

**Challenge:** The base class `PsiThermo::correct()` was being called instead of the custom two-temperature version.

**Solution:** Modified ShockThermo with pointer to custom thermodynamics:
```cpp
highEnthalpyMulticomponentThermo* heThermoPtr_;

// Inside thermophysicalPredictor():
if (customThermoActive)
    heThermoPtr_->correct_he();
else
    thermo_.correct();
```
This ensures correct two-temperature energy update at each time step.

### Mutation++ Integration

**Thermochemical calculations via Mutation++:**
- Relaxation time scales
- Chemical kinetics integration
- Thermodynamic properties

**Critical requirement:** Species ordering in OpenFOAM must match Mutation++ library order.

**Core call:**
```cpp
mutationMixPtr_->step(dt, rho, Y, Et, Ev, Ttr, Tv);
```
This updates both T and Tᵥ consistently with energy modes.

## Results

### 1. Temperature Relaxation (Single Cell, Zero-D)

**Initial conditions:**
- Tt = 12000 K (translational)
- Tv = 2000 K (vibrational)

**Observed behavior:**
- Tt **decreases rapidly** (energy transfer to vibrational mode)
- Tv **increases** (receives energy from translation)
- Both **converge to equilibrium temperature**

✅ **Validates correct vibration–translation relaxation physics**

The temperature intersection point corresponds to energy mode balancing.

### 2. Multi-Cell Speedup Analysis

**Domain:** 10 × 10 × 10 grid (1000 cells)

**Parallelization strategy:** MPI

**Key findings:**
- Significant speedup observed up to **4 processes**
- Performance plateaus after 4 processes (MPI overhead dominates)
- Hardware limitation: 4 performance cores on test machine
- OpenMP: Not beneficial (single-core improvement only)

**Conclusion:** MPI scaling works well and is recommended for larger simulations.

### 3. Validation with Reference Paper

Casseau et al. (2016) comparison:

| Parameter | Result | Status |
|-----------|--------|--------|
| Tt evolution | Excellent agreement | ✅ |
| Tv evolution | Excellent agreement | ✅ |
| Species evolution | Correct reproduction | ✅ |
| O₂ consumption | Accurate | ✅ |
| O production | Accurate | ✅ |
| Chemical kinetics | Properly integrated | ✅ |

**Overall:** Model **successfully validated** against literature.

## Model Validation

The implementation is validated in three aspects:

1. **Physics correctness** – Temperature relaxation reproduces expected behavior
2. **Numerical accuracy** – Multi-cell simulations show proper MPI scalability  
3. **Literature comparison** – Results agree with Casseau et al. (2016)

The framework is suitable for **high-enthalpy hypersonic reacting flow simulations**.

## Contributing

Contributions are welcome! Please ensure:
- Code follows OpenFOAM conventions
- New features include validation tests
- Documentation is updated
- MPI parallelization is verified for scalability

## References

**Primary Reference:**
> Casseau, V., et al. (2016). "A two-temperature open-source CFD model for hypersonic reacting flows." Computer Physics Communications, 203, 104-120.

**Related Work:**
- OpenFOAM Foundation documentation
- Mutation++ library documentation
- Hypersonic aerothermodynamics literature

## License

[Specify your license here - e.g., GPL-3.0, MIT, etc.]

## Authors & Contributors

Project developed as part of advanced CFD coursework.

Contributors:
- Allanda
- Bergamaschi
- Esposito
- Faggion
- Grassi
- Venezia
- Virgulti

## Contact & Support

For questions or issues:
1. Check the full project report: `ReportMod2_*.pdf`
2. Review tutorials in `tutorials/` directory
3. Consult Mutation++ documentation
4. Open an issue on GitHub

---

**Last Updated:** March 2026  
**Status:** Validated and tested ✅