# Spectral Semi-Lagrangian Semi-Implicit Shallow-Water Model Research Code

This repository contains Fortran research code for the global shallow-water equations 
on the sphere using a spectrally semi-Lagrangian semi-implicit (SISL) formulation. 
The code supports both double Fourier sphere (DFS) and spherical harmonic transform (SHT) methods for
computing the spatial derivatives, as implemented by [Yoshimura](https://doi.org/10.5194/gmd-15-2561-2022-supplement). 
The main novelty of this model is that it uses spectrally accurate interpolation for the semi-Lagrangian part via DFS, which differs from standard SISL models that use low-order interpolation. To make the DFS interpolation compuationally competitive, we use a non-uniform fast Fourier transform (NUFFT), as implemented in [FINUFFT](https://github.com/flatironinstitute/finufft).

The reposistory includes both DFS and SHT formulations, local transform-library dependencies, build scripts, and
batch-run scripts for several standard test cases in the literature, including Williamson et. al. test cases 1, 2, 5, & 6, and the Galewsky et. al. test case.

The main model code is split into two formulations:

- `shallow_water_SH/`: spherical-harmonic transform version.
- `shallow_water_DFS/`: double-Fourier-series version.

The directories `bihar/` and `ispack-3.0.1/` provide the supporting
transform libraries used by the model builds.

The repository modifies the supplementary code associated with [Yoshimura
(2022)](https://doi.org/10.5194/gmd-15-2561-2022-supplement).

The bundled support libraries should also be cited when they are used:
`bihar/` contains FFTPACK-style double-precision real, sine, and cosine
transform routines, for which the standard reference is Swarztrauber
(1982). `ispack-3.0.1/` is ISPACK3 by Keiichi Ishioka; see the bundled
`ispack-3.0.1/README` and `ispack-3.0.1/COPYRIGHT`.

## Directory Layout

| Path | Purpose |
| --- | --- |
| `shallow_water_SH/` | Spherical-harmonic shallow-water and advection models, transform utilities, diagnostics, and tests. |
| `shallow_water_DFS/` | Double-Fourier-series shallow-water and advection models, DFS transforms, derivative operators, elliptic solvers, diagnostics, and tests. |
| `bihar/` | Static library of FFTPACK-style real sine/cosine/FFT routines. Builds `libbihar.a`. |
| `ispack-3.0.1/` | Bundled ISPACK3 spectral-transform library. Builds `libispack3.a`; see its own `README` and package docs. |
| `docs/CODE_STRUCTURE.md` | Developer-oriented map of entry points, modules, outputs, and common configuration constants. |

## External Dependencies

The Makefiles assume the following tools and libraries are available:

- `gfortran`
- OpenMP support for the Fortran compiler
- FFTW3 and FFTW3 OpenMP libraries: `libfftw3`, `libfftw3_omp`
- BLAS/LAPACK for the spherical-harmonic build
- FINUFFT installed under `$(HOME)/finufft`
- A C/C++ runtime suitable for the linked FINUFFT and FFTW libraries

The current Makefiles contain machine-specific paths such as
`$(HOME)/finufft/include`, `$(HOME)/finufft/lib`, and
`$(HOME)/finufft/lib-static`. Adjust these paths if FINUFFT is installed
elsewhere.

## Build Order

Build the support libraries before building either shallow-water model:

```sh
cd bihar
make

cd ../ispack-3.0.1
make
```

Then build one of the model directories:

```sh
cd ../shallow_water_SH
make
```

or:

```sh
cd ../shallow_water_DFS
make
```

Use `make clean` in each directory to remove compiled objects, modules,
and executables created by that directory.

## Executables

`shallow_water_SH/` builds:

- `sw_sh`: shallow-water model using spherical harmonics.
- `adv_sh`: passive/advection test case using spherical harmonics.
- `test_laplacian_helmholtz`: operator accuracy test for Laplacian and
  Helmholtz inversion.
- `test_truncate`: transform/truncation test.

`shallow_water_DFS/` builds:

- `sw_dfs`: shallow-water model using double Fourier series.
- `adv_dfs`: passive/advection test case using double Fourier series.
- `test`: DFS transform/operator test driver.

## Running

The forecast scripts set OpenMP runtime variables and then launch the
corresponding executable:

```sh
cd shallow_water_SH
./fcst_sw_sh.sh
./fcst_adv_sh.sh

cd ../shallow_water_DFS
./fcst_sw_dfs.sh
./fcst_adv_dfs.sh
```

The scripts also include PJM scheduler directives for batch execution on
systems that support PJM. On a local workstation, the `#PJM` lines are
comments and the scripts can still be run directly after the executables
are built.

For local runs, reduce `OMP_NUM_THREADS` and `OMP_STACKSIZE` if the machine
does not have the memory assumed by the batch settings.

## Configuration

Most experiments are configured by editing `integer, parameter` and
`real(8), parameter` constants near the top of each program file:

- `JCN_INITIAL`: selects the initial condition or benchmark case.
- `JCN_HDIFF`: selects no diffusion, second-order diffusion, or fourth-order
  hyperdiffusion where supported.
- `JCN_MONIT`: enables monitor output.
- `JCN_MONIT_SPECTRUM`: enables kinetic-energy spectrum output where present.
- `INTHR_MONIT`: monitor interval in model hours.
- `IT`: interpolation method for semi-Lagrangian departure-point fields.
  Use `IT = 1` for spectral interpolation and `IT = 0` for low-order
  Lagrange interpolation.
- `IMAX`, `JMAX`, `NMAX`: horizontal grid and spectral truncation.
- `TIMESTEP`: model time step in seconds.

DFS programs also expose:

- `JCN_DFS`: DFS basis/coefficient method.
- `JCN_GRID`: latitude grid choice.
- `JCN_ZONALFILTER` and `M0`: zonal Fourier filter settings.

After changing these constants, rebuild the affected executable.

## Output Files

Diagnostics are written as direct-access binary `.dr` files with matching
GrADS control files (`.ctl`). Common outputs include:

- `data.dr`: gridded monitor fields.
- `norm.dr`: norm/error diagnostics.
- `spectrum.dr`: kinetic-energy spectrum diagnostics when enabled.
- `weight_lat.dr`: latitudinal quadrature weights.

The `grads` modules write control files describing the grid, levels, and
time axis for visualization tools that understand GrADS control files.

## Converting Output To NetCDF

The model writes GrADS datasets as a pair of files: a binary data file
(`.dr`) and a control file (`.ctl`). Convert these datasets to NetCDF with
Climate Data Operators (CDO):

```sh
cdo -f nc import_binary data.ctl data.nc
```

Run the command from the directory containing both `data.ctl` and
`data.dr`. The same pattern applies to other outputs:

```sh
cdo -f nc import_binary norm.ctl norm.nc
cdo -f nc import_binary spectrum.ctl spectrum.nc
cdo -f nc import_binary weight_lat.ctl weight_lat.nc
```

To convert every GrADS control file in a run directory:

```sh
for ctl in *.ctl; do
  base=${ctl%.ctl}
  cdo -f nc import_binary "$ctl" "$base.nc"
done
```

If CDO reports a missing data file, check the `DSET` line inside the `.ctl`
file. It should point to the matching `.dr` file, usually with a relative
path such as `DSET ^data.dr`.

## Notes For Maintainers

- This workspace is not arranged as one monolithic build system. Each
  directory has its own Makefile and cleanup target.
- The shallow-water programs are intentionally parameter-driven Fortran
  drivers rather than command-line applications. Experiment changes usually
  require editing constants and recompiling.
- `ispack-3.0.1/` is third-party bundled source. Prefer documenting local
  integration points rather than editing ISPACK internals unless you are
  intentionally modifying the library.

## References

- Yoshimura, H.: Improved double Fourier series on a sphere and its
  application to a semi-implicit semi-Lagrangian shallow-water model,
  Geosci. Model Dev., 15, 2561-2597, 2022. [Paper](https://doi.org/10.5194/gmd-15-2561-2022) &
  [Supplement](https://doi.org/10.5194/gmd-15-2561-2022-supplement).
- D. L. Williamson, J. B. Drake, J. J. Hack, R. Jakob, and P. N. Swarztrauber:
  A standard test set for numerical approximations to the shallow water equations
  in spherical geometry, J. of Comp. Phys., 102, 1992.
- J. Galewsky, R. K. Scott, and L. M. Polvani: An initial-value problem for
  testing numerical models of the global shallow-water equations, Tellus A: Dyn.
  Meteorol. Oceanogr., 56, 429–440, 2004.
- A. H. Barnett, J. Magland, and L. af Klinteberg: A parallel nonuniform fast
  Fourier transform library based on an “exponential of semicircle” kernel,
  SIAM J. Sci. Comput., 41, C479–C504, (2019).
- A. H. Barnett, J. Magland, and L. af Klinteberg: Flatiron Institute nonuniform
  fast Fourier transform libraries (FINUFFT). [code](https://github.com/flatironinstitute/finufft).
- Swarztrauber, P. N.: Vectorizing the FFTs, in: Parallel Computations,
  edited by G. Rodrigue, Academic Press, 51-83, 1982. This is the
  reference listed by Netlib for FFTPACK, the FFT package whose real,
  sine, and cosine transform routines correspond to the support routines
  bundled in `bihar/`.
- Ishioka, K.: ISPACK library, GFD-DENNOU Club [code](https://www.gfd-dennou.org/arch/ispack/),
  cited by the bundled ISPACK README for `ispack-3.0.1/`.
- Ishioka, K.: A new recurrence formula for efficient computation of
  spherical harmonic transform, J. Meteorol. Soc. Jpn., 96, 241-249, 2018.
  [Paper](https://doi.org/10.2151/jmsj.2018-019).
