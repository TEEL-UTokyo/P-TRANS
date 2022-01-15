# P-TRANS
![logo](http://www.phonon.t.u-tokyo.ac.jp/p-trans-jp/_images/PTRANS_logo.png)
## Overview
> `P-TRANS` is a software to simulate phonon transport in arbitrary nanostructures based on the Boltzmann transport theory. 
> The core part of the software is a Monte Carlo ray-tracing solver that uses the stochastic method to sample phonon transport in various nanostructures. 
> The solver is written in FORTRAN with well-organized modules utilizing the OpenMP and MPI parallelization to achieve high performance. 
> The software can handle a typical calculation of nanostructures with a feature size of ∼1000 nm on a normal PC within seconds, which is useful when exploring the wide variety of structures, for instance for structural optimization or high-throughput screening in material informatics.
> The software takes the bulk phonon properties obtained from the first-principles anharmonic lattice dynamics calculations as inputs. 
> The database of bulk phonon properties is available on the software webpage, which currently contains ∼60 materials and will increase in the future. We also provide a graphical user interface to design the nanostructures and to run the simulation for non-expert users. In this paper, we introduce the theoretical background, the technical details and implementation, and the applications of the software. 

## Homepage
[P-TRANS](http://www.phonon.t.u-tokyo.ac.jp/p-trans/)

## Reference
- C. Shao, T. Hori and J. Shiomi*, “P-TRANS: A Monte Carlo ray-tracing solver for phonon transport in nanostructures”, Computer Physics Communications, under review, minor revision. 
