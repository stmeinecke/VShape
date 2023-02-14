# VShape
C++ code to simulate a V-shaped mode-locked semiconductor disk laser

## What

The code numerically integrates a minimal delay-differential equation model for a V-shaped mode-locked semiconductor disk laser via a traditional fourth-order Runge-Kutta solver with cubic Hermite spline interpolation for the delay array. A number of heuristics are are provided to characterize the resulting time series.

## Setup
The code can be compiled via the make tool and the (minimal) makefile. The utilized fftw library is provided in the /incl folder.

The code has been tested on a Intel i7-8700 running Debian Bullseye (compiled with gcc-10 (10.2.1-6)).

## Usage
The simulation code is controlled via command line inputs. The relevant flags can be found in the VShape.cpp file (main file).

For example, the input
```
./VShape -TS -intTime 20000 -nout 10
```
simulates a time series (`-TS`) with the default system parameters and initial conditions for 20000 time units (`-intTime 20000`) where only every tenth time step is saved (`-nout 10`). Note that only the last 100 time units are saved (can be controlled via `-outTime`). Once finished, the resulting time trace is saved to 'data/out_TS' and a number of time series characteristics are computed and printed to the terminal.

## Notes

The simulation code represents a typical small-scale scientific Frankenstein project. Many parts have been reused and repurposed from other projects (as indicated with comments in the respective files). Hence, some sections and features may appear confusing and unclear at first.

In particular, the time-series evaluation class evalMLTS is a monster that has grown for a couple of years. It collects a number of different heuristics that have been developed for different mode-locked laser models. Most of them require thresholds to be set to reasonable values. This is the users responsibility.

The delay-differential equations (DDE) solver, on the other hand, has been thoroughly tested and proven to be reliable and robust for a number of different projects. Note, that the traditional fourth-order Runge-Kutta method requires the user to set an appropriate time step and ensure sufficient convergence.

## References

This code (and prior versions) has been used for the following publications:

> Meinecke, S. and Lüdge, K., Optimizing the Cavity-Arm Ratio of V-Shaped Semiconductor Disk Lasers, Phys. Rev. Appl. 18, 064070, (2022). [[doi](https://doi.org/10.1103/PhysRevApplied.18.064070)]

> Meinecke, S., Spatio-Temporal Modeling and Device Optimization of Passively Mode-Locked Semiconductor Lasers, Springer Theses, Springer, Cham, ISBN 978-3-030-96247-0, (2022). [[doi](http://dx.doi.org/https://doi.org/10.1007/978-3-030-96248-7)]

> Hausen, J., Meinecke, S., Javaloyes, J., Gurevich, S. V. and Lüdge, K., Phase-Incoherent Photonic Molecules in V-Shaped Mode-Locked Vertical-External-Cavity Surface-Emitting Semiconductor Lasers, Phys. Rev. Appl. 14, 044059, (2020). [[doi](http://dx.doi.org/https://doi.org/10.1103/physrevapplied.14.044059)]
 
> Hausen, J., Meinecke, S., Lingnau, B. and Lüdge, K., Pulse Cluster Dynamics in Passively Mode-Locked Semiconductor Vertical-External-Cavity Surface-Emitting Lasers, Phys. Rev. Appl. 11, 044055 (2019). [[doi](http://dx.doi.org/10.1103/physrevapplied.11.044055)]

> Hausen, J., Meinecke, S. and Lüdge, K., Bifurcation scenario leading to multiple pulse emission in VECSELs, Proc. SPIE 10901, 109010F (2019). [[doi](http://dx.doi.org/10.1117/12.2513751)]

> Waldburger, D., Alfieri, C. G. E., Link, S. M. , Meinecke, S., Jaurigue, L. C., Lüdge, K. and Keller, U., Multipulse instabilities of a femtosecond SESAM-modelocked VECSEL, Opt. Express 26, 17, 21872 (2018). [[doi](http://dx.doi.org/10.1364/oe.26.021872)]
