# Bayesian-Rating-Curve
This repository contains the Bayesian rating model codes used to analyze rating curve dynamics and heteroscedastic flow error modeling for five hydrometric stations in India.
The codes are developed based on the fundamentals provided by the developers of the BaRatin Framework (Le Coz et al., 2014).

## Prerequisites
To run these models, you must have R installed along with the 'rstan' package and the associated dependencies. 

### Installing rstan
The `rstan` package is the R interface to Stan. For comprehensive installation instructions and platform-specific troubleshooting, please refer to the official [RStan Getting Started guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

```r
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
```
#### Reference
Le Coz, J., B. Renard, L. Bonnifait, F. Branger, and R. Le Boursicaud. (2014). “Combining hydraulic knowledge and uncertain gaugings in the estimation of hydrometric rating curves: A Bayesian approach.” Journal of Hydrology, 509: 573–587. https://doi.org/10.1016/j.jhydrol.2013.11.016.
