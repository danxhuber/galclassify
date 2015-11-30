# galclassify
Perform stellar classifications using a galactic population synthesis model

###Documentation:
- Sharma et al. (2011): http://adsabs.harvard.edu/abs/2011ApJ...730....3S
- Huber et al. (2016, in prep)


###Dependencies:
* Galaxia: http://galaxia.sourceforge.net/ (the galaxia executable needs to be in your path) <br/>
* IDL/GDL: ASTROLIB, Coyote Library <br/> 
* Code is tested and functional with GDL 0.9.5 (open source implementation of IDL; using GDL additionally requires installation of the CMSVLIB library: https://www.physics.wisc.edu/~craigm/idl/cmsave.html to restore IDL save files)


###classifyk2:

Contains code used in Huber et al. (2016) to generate stellar properties posteriors for ~120,000 stellar sources targeted by the K2 mission. 

[under construction]


###galclassify:

Stand-alone version to classify any given star (not necessarily observed by K2) given some input observables. Unlike classifyk2 this code generates a synthetic population for each star given some radius, rather than using a precomputed population over a large field. Suitable for individual stellar classifications, rather slow if you want to classify large samples of stars.

[under construction]
