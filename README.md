# CoExpansionValidation
An R package for validating ABC co-expansion inference 



Tests of synchronous community expansion are commonly conducted using hierarchical Approximate Bayesian Computation (hABC), a statistical framework for inferring the degree of concordance across species. However, this framework is often used without demonstrating adequate performance.

The examplerun file in the vignettes will go through step-by-step an example of how to generate one pseudo-observed dataset, and use ABC to "infer" the number of co-expansion events. Obviously, to apply these performance assessment in a empirical study, many replicates of pseudo-observed datasets are needed.   

To install this package through use devtools:

```r
devtools::install_github("huatengh/CoExpansionValidation", upgrade_dependencies = TRUE)
```
