
# Finding Eigenpairs of Interest Using Ginkgo

GPU code for finding eigenpairs of large sparse matrices near randomly chosen energies, without full diagonalization, using Ginkgo. The code solves for the eigenpairs and reports their accuracy; computing expectation values or other quantities from the eigenvectors is left to the user. This repository accompanies the paper:

> N. A. Chlanda, J. Keim, I. Martinez-Robles, S. E. Spielman, T. J. Carroll, and M. W. Noel, *A random sampling method to test the weak Eigenstate Thermalization Hypothesis* (2026). [[arXiv link / DOI](https://arxiv.org/abs/2610.30070)]

## Table of Contents
1. [The Method](#the-method)
2. [Repository Structure](#repository-structure)
3. [Requirements](#requirements)
4. [Building and Running](#building-and-running)
5. [Output](#output)

## The Method

Given a sparse real symmetric matrix $H$, the code:

1. Estimates the spectral bounds of $H$ with a short Lanczos run and rescales $H$ onto $[-1, 1]$.
2. Estimates the density of states $\rho(E)$ with the kernel polynomial method (KPM).
3. Draws energies at random from a central fraction of the spectrum, using by $\rho(E)$.
4. Opens a window around each energy expected to contain a set number of eigenvalues.
5. Finds the eigenpairs in that window with the main loop: Chebyshev window filter → SVQB orthonormalization → Rayleigh–Ritz.
6. Keeps the in-window Ritz pairs with relative residual below $10^{-3}$ and repeats until the requested number of eigenpairs has been collected.

The KPM and CFD routines can also be used on their own.

## Repository Structure

```
├── FindingEigenpairsOfInterest.cpp   # main source: Lanczos, KPM, CFD, and sampling loop
├── Makefile                          # build file
└── README.md
```

## Requirements

- An NVIDIA GPU with CUDA
- [Ginkgo](https://github.com/ginkgo-project/ginkgo) built with CUDA support
- Intel MKL (LAPACKE)
- NetCDF (C library)

Edit the paths in the `Makefile` to match your installation.

## Building and Running

### Build
```bash
git clone https://github.com/NickChlanda/Finding-Eigenpairs-of-Interest-Using-Ginkgo.git
cd Finding-Eigenpairs-of-Interest-Using-Ginkgo
make
```

### Run
```bash
./FindingEigenpairsOfInterest <matrix.mtx> <gpu_id> <NP> <NT> <max_outer> <target_frac> <mc_target> <ns_max>
```

| Argument      | Description                                                        |
|---------------|--------------------------------------------------------------------|
| `matrix.mtx`  | Sparse matrix in MatrixMarket coordinate format                    |
| `gpu_id`      | GPU device to run on                                               |
| `NP`          | Chebyshev polynomial filter degree                                 |
| `NT`          | Target number of eigenvalues per window                            |
| `max_outer`   | Maximum CFD iterations for reach window                            |
| `target_frac` | Fraction of the spectrum to sample from (e.g. `0.05` = middle 5%)  |
| `mc_target`   | Total number of eigenpairs to collect                              |
| `ns_max`      | Maximum number of search vectors (the code uses `min(2*NT, ns_max)`)|

Example:
```bash
./FindingEigenpairsOfInterest H_N13.mtx 0 3500 150 8 0.05 2000 300
```

The parameters used in the paper are listed in Table I. For similar sized Hamiltonians to ours, I'd start with our parameters and change as needed.

## Output

- **`spectrum.txt`**: the KPM density of states (energy, normalized $\rho$).
- **`ritz_hits_gpu<gpu_id>.nc`**: one NetCDF variable `DATA` that stores three values per accepted eigenpair: the eigenvalue, the residual norm, and the first component of the eigenvector. The full eigenvector is held in memory (`ap.vec`) and can be written out if needed.
- **Console**: progress, per-iteration residuals, and a final summary of the eigenpair accuracy (mean, best, and worst relative residual), duplicates removed, and run time.
