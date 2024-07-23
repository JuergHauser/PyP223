# PyP223

The P223F software suite (Raiche et al.2007) is an extensive collection of tools for EM forward modelling an inversion, that in the past has been made availablee here https://p223suite.sourceforge.net/

This repository contains ctypes based Python wrappers for the variant of LeroiAir employed in Hauser et al. 2016, specfically the relevant subset of the functions in libp223 that is built as part of wiglaf and used solve the forward problem. Thus the files, in pyp223/p223 are overwritten whenever there is an update to https://gitlab.com/jrh/wiglaf. 

The version of LeroiAir in wiglaf is subtly different from previously released versions:
- The code-base has been re-modularised to avoid name clashes in libp223 with other P223F subroutines/programs
- There is a functional interface to forward model specific systems
- Plates are parametrised differently, see Hauser et. al (2016) 

Currently the Python API provides two functions to model the response of a Tempest and VTemMax system to a thin plate in the halfspace of a layered
earth

```
formod_tempest_data(nlyr,nstat,...)
formod_vtem_max_data(nlyr,nstat,...)
```

See `demos/verification` for the comparisons between  results obtained with the functions provided here and computations perforemd with LeroiAir.

## Installation
```
pip install .
```

## Licensing
PyP223 is released as  GPL-3.0-or-later.

### References
Raiche, A., Sugeng, F. and Wilson, G. (2007) Practical 3D EM inversion the P223F software suite, ASEG Extended Abstracts, 2007:1, 1-5

Hauser, J., Gunning, J., & Annetts, D. (2016). Probabilistic inversion of airborne electromagnetic data for basement conductors. Geophysics, 81(5), E389-E400.
