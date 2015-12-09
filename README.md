# galclassify
Perform stellar classifications using a synthetic stellar population generated by Galaxia

###Documentation:
- Sharma et al. (2011): http://adsabs.harvard.edu/abs/2011ApJ...730....3S
- Huber et al. (2015, submitted)


###Dependencies:
* Galaxia: http://galaxia.sourceforge.net/ (the galaxia executable needs to be in your path) <br/>
* IDL/GDL: ASTROLIB, Coyote Library <br/> 
* Code is tested and functional with GDL 0.9.5 (open source implementation of IDL; using GDL additionally requires installation of the CMSVLIB library: https://www.physics.wisc.edu/~craigm/idl/cmsave.html to restore IDL save files)


##galclassifyk2:

Code used to classify stellar sources targeted by the K2 mission, which can be used to e.g. extract posteriors for any given EPIC target or derive classifications for a new campaign. Code will be added after peer review of the paper (submitted to ApJS).


##galclassify:

Stand-alone version of the Galaxia classification code to classify any given star (not necessarily observed by K2) given some input observables. Unlike classifyk2 this code generates a synthetic population for each star given some radius, rather than using a precomputed population over a large field. Suitable for individual stellar classifications, but slow if you want to classify large samples of stars.

####Usage:

Calling sequence from IDL/GDL:
```
galclassify,input=input
```

**Required input parameters:**	<br/>
input ... plain ascii file containing observables in a fixed format; see tychostar.txt for an example <br/>

**Optional input parameters:**	<br/>
sample	... specifies how many posterior samples should be saved in a text file  <br/>
pl      ... set if you want to see plots

**Output:**	<br/>
If sample is not set, the default output will be the median and 1-sigma confidence interval for each parameter (printed in the terminal) and graphical output of posteriors and observables. If sample is set, the text files containing the posteriors (both discrete and sampled) will be written to the output/ directory

**Example:** <br/>
Classify some random Tycho star:
```	
galclassify,input="tychostar.txt"
```

Classify some random Tycho star and save posteriors with 10000 samples:
```	
galclassify,input="tychostar.txt",sample=10000
```

Classify the Sun at J=10mag and plot result:
```	
galclassify,input="sun_at_J10.txt",/pl
```

####Random and important notes:
* For bright stars Galaxia will take a very long time to run, and most likely the population will still be too sparse for meaningful inference - in this case it's better to use traditional isochrone modeling.
* The posterior sampling is still very crude, and for stars with poor model coverage the samples will be dominated by a few distinct models
