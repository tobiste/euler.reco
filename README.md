<!-- badges: start -->
[![R-CMD-check](https://github.com/tobiste/euler_reco/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tobiste/euler_reco/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# euler.reco
Algorithms to find and evaluate the Euler pole solution describing the 
orientation of geological structures

## Method
The small circle of best fit and the corresponding Euler pole
(axis) of rotation for an array of points are established by treating each 
point as the locus of a unit vector from the center of Earth and by
minimizing the sums of the squares of the deviations of these unit vectors from
a circular cone surface with its apex at the center of Earth 
(Ramsay, 1967, p. 20-21).

## Prerequisites

You must have R installed on your system (see http://r-project.org). To install **euler.reco** from Github, you also need the `remotes` package. This can be installed by typing the following code at the R command line prompt:

```
install.packages("remotes")
```

## Installation

The most recent development version of `euler.reco` is available from Github and can be installed on your system as follows:

```
remotes::install_github('tobiste/euler.reco')
library('euler.reco')
```

## Documentation

The detailed documentation can be found at 
https://tobiste.github.io/euler.reco/


## Author
Tobias Stephan

## How to cite
...

## License
MIT License

## References
- <div class="csl-entry">Ramsay, J. G. (1967). <i>Folding and fracturing of rocks. McGraw-Hill Book Company</i>, New York</div>

- <div class="csl-entry">Gray, N.H., Geiser, P.A., Geiser, J.R. (1980). On the least-square fit of small and great circles to spherically projected data. <i>Mathematical Geology</i>, <i>12</i>(3)</div>

- <div class="csl-entry">Kroner, U., Roscher, M., Romer, R. L. (2016). Ancient plate kinematics derived from the deformation pattern of continental crust: Paleo- and Neo-Tethys opening coeval with prolonged Gondwana–Laurussia convergence. <i>Tectonophysics</i>. <i>681</i>, 330-233. https://doi.org/10.1016/j.tecto.2016.03.034</div>


- <div class="csl-entry">Wdowinski, S. (1998). A theory of intraplate tectonics. <i>Journal of Geophysical Research: Solid Earth</i>, <i>103</i>(3), 5037–5059. http://dx.doi.org/10.1029/97JB03390</div>

- <div class="csl-entry">Price, R. & Carmicheal, D. (1986). Geometric test for Late Cretaceous-Paleogene intracontinental transform faulting in the Canadian Cordillera. <i>Geology</i>, 468, 14(6), 468-471. doi: 10.1130/0091-7613(1986)14<468:GTFLCI>2.0.CO;2</div>
