# GOTHiC++

GOTHiC++ is a C++ implementation of the Bioconductor R package [GOTHiC](https://doi.org/doi:10.18129/B9.bioc.GOTHiC "GOTHiC on Bioconductor")

This is Hi-C analysis software which uses a cumulative binomial test to detect interactions between distal genomic loci that have significantly more reads than expected by chance in Hi-C experiments. It takes mapped paired NGS reads as input and gives back the list of significant interactions for a given bin size in the genome.

**GOTHiC Authors**: Borbala Mifsud and Robert Sugar
**GOTHiC Citation**: Mifsud B, Sugar R (2020). GOTHiC: Binomial test for Hi-C data analysis. R package version 1.24.0.

**GITHiC++ Authors**: Richard Thompson and Borbala Mifsud
**Contact Email**: ithompson[at]hbku.edu.qa

###Requirements
+ Cmake 3.4
+ C++17-compatible compiler

###Building and Installing

In order to compile GOTHiC++, the authors advise creating a build folder within the root:

```bash
mkdir build && cd build
```

Compiling GOTHiC++ follows a typical cmake pipeline:

```bash
cmake ../source
make
```


###Testing

_To Code!_
 