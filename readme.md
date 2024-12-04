# PyP223

The P223F software suite (Raiche et al.2007) is an extensive collection of tools for EM forward modelling and inversion, that in the past has been made available on a website maintained by AMIRA. While the original site no longer exists, it and the related files have been captured by the wayback machine and are available [here](https://web.archive.org/web/20160313045828/http://amirainternational.com/web/site.asp?page=projectpages/p223f_software&section=news)

Several repositories exist that provide variants of P223F, seeking to re-enable features originally disabled by Art Raiche, providing a build system etc...

[https://github.com/dwannetts/P223_Public](https://github.com/dwannetts/P223_Public)

[https://github.com/KoBoldMetals/P223_Public](https://github.com/KoBoldMetals/P223_Public)

They all appear to adhere to the orignal philosophy, that is read one or several input files from disk and then write the results to one or several output files.

This repository contains ctypes based Python wrappers for the variant of LeroiAir employed in Hauser et al. (2016), specfically the relevant subset of the functions/files from libp223 that is built as part of wiglaf and used solve to the forward problems. Thus the files, in pyp223/p223 are currently overwritten whenever there is an update to [wiglaf](https://gitlab.com/jrh/wiglaf). 

The version of LeroiAir in libp223 that is part of wiglaf is subtly different from the Version provided by P223F.
- The code-base has been re-modularised to avoid name clashes in libp223 which includes other P223F functionality from Airbeo, Beowulf and Leroi.
- There is a functional interface to forward model specific systems
- Plates are parametrised differently, see Hauser et. al (2016) 

Currently the Python API provides two functions to model the response of a TEMPEST and VTEM max system to a thin plate in the halfspace of a layered
earth. It would be straightforward to expose the functions in `p223/clibleroiair.f90` as well as `clibairbeo.f90`, `clibbeowulf.f90`  and `clibleroi.f90` from wiglaf's libp223 for other airborne and ground systems by adding the necessary python code in `/src/pyp223/_pyp223.py`. 

```
formod_tempest_data(nlyr,nstat,...)
formod_vtem_max_data(nlyr,nstat,...)
```

See `demos/verification` for the comparisons between results obtained with the functions provided here and computations performed with LeroiAir from P223F. 

An example for how the Python wrappers can be used for the inversion of VTEM max data using [CoFI](https://inlab.au/cofi/) is given [here](https://github.com/inlab-geo/cofi-examples/blob/main/tutorials/thin_plate_inversion/thin_plate_inversion.ipynb).


## Installation
```
git clone https://github.com/JuergHauser/PyP223.git
cd PyP223
pip install .
```

## Licensing
PyP223 is released as  GPL-3.0-or-later.

### References
Raiche, A., Sugeng, F. and Wilson, G. (2007) Practical 3D EM inversion the P223F software suite, ASEG Extended Abstracts, 2007:1, 1-5

Hauser, J., Gunning, J., & Annetts, D. (2016). Probabilistic inversion of airborne electromagnetic data for basement conductors. Geophysics, 81(5), E389-E400.
