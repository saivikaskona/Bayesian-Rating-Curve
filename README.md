# Bayesian-Rating-Curve
This repository contains the Bayesian rating model codes used to analyze rating curve dynamics and heteroscedastic flow error modeling for five hydrometric stations in India.
The codes are developed based on the fundamentals provided by the developers of the BaRatin Framework (Le Coz et al., 2014).

## Prerequisites
To run these models, you must have R installed along with the 'rstan' package and the associated dependencies. 

### Installing rstan
The `rstan` package is the R interface to Stan. You can install it by running the following command in your R console:
```r
install.packages("rstan", repos = "[https://cloud.r-project.org/](https://cloud.r-project.org/)", dependencies = TRUE)


## References
If you use these codes or the BaRatin framework in your research, please cite the following original work:

* **Le Coz, J., Renard, B., Bonnifait, L., Branger, F., & Le Boursicaud, R. (2014).** *Combining hydraulic knowledge and uncertain gaugings in the estimation of hydrometric rating curves: A Bayesian approach.* Journal of Hydrology, 509, 573-587. [https://doi.org/10.1016/j.jhydrol.2013.11.016](https://doi.org/10.1016/j.jhydrol.2013.11.016)
