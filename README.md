# VShape
C++ code to simulate a V-shaped mode-locked semiconductor disk laser

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



## References

This code (and prior versions) has been used for the following publications:

> Meinecke, S., Spatio-Temporal Modeling and Device Optimization of Passively Mode-Locked Semiconductor Lasers, Springer Theses, Springer, Cham, ISBN 978-3-030-96247-0, (2022). [[doi](http://dx.doi.org/https://doi.org/10.1007/978-3-030-96248-7)]

> Hausen, J., Meinecke, S., Javaloyes, J., Gurevich, S. V. and Lüdge, K., Phase-Incoherent Photonic Molecules in V-Shaped Mode-Locked Vertical-External-Cavity Surface-Emitting Semiconductor Lasers, Phys. Rev. Appl. 14, 044059, (2020). [[doi](http://dx.doi.org/https://doi.org/10.1103/physrevapplied.14.044059)]
 
> Hausen, J., Meinecke, S., Lingnau, B. and Lüdge, K., Pulse Cluster Dynamics in Passively Mode-Locked Semiconductor Vertical-External-Cavity Surface-Emitting Lasers, Phys. Rev. Appl. 11, 044055 (2019). [[doi](http://dx.doi.org/10.1103/physrevapplied.11.044055)]

> Hausen, J., Meinecke, S. and Lüdge, K., Bifurcation scenario leading to multiple pulse emission in VECSELs, Proc. SPIE 10901, 109010F (2019). [[doi](http://dx.doi.org/10.1117/12.2513751)]

> Waldburger, D., Alfieri, C. G. E., Link, S. M. , Meinecke, S., Jaurigue, L. C., Lüdge, K. and Keller, U., Multipulse instabilities of a femtosecond SESAM-modelocked VECSEL, Opt. Express 26, 17, 21872 (2018). [[doi](http://dx.doi.org/10.1364/oe.26.021872)]
