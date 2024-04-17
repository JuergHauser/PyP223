# Pyp223

Ctypes based python wrappers for the variant of LeroiAir employed in Hauser et al. 2016, specfically the relevant subset of the functions in libp223 that is built as part of wiglaf and used solve the forward problem. Thus the files, in pyp223/modules are overwritten whenever there is an update to https://gitlab.com/jrh/wiglaf.

The version of LeroiAir in wiglaf is subtly different from previously released versions:
- The code-base has been re-modularised to avoid name clashes with other P223 
  subroutines/programs
- There is a functional interface to forward model specific systems
- Plates are parametrised differently, see Hauser et. al (2016) 
- There is an ability to model plunge of a thin plate

Currently the python API provides two functions to model the response of a Tempest and VTemMax system to a thin plate in the halfspace of a layered
earth

```
formod_tempest_data(nlyr,nstat,...)
formod_vtem_max_data(nlyr,nstat,...)
```

## Installation
```
pip install .
```

## Licensing
PyP223 is released under the GNU General Public License version 2.

### References
Hauser, J., Gunning, J., & Annetts, D. (2016). Probabilistic inversion of airborne electromagnetic data for basement conductors. Geophysics, 81(5), E389-E400.
