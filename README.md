# R package FindIt [![Build Status](https://travis-ci.org/kosukeimai/FindIt.svg?branch=master)](https://travis-ci.org/kosukeimai/FindIt)  [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/FindIt)](https://cran.r-project.org/package=FindIt)

## FindIt: Finding Heterogeneous Treatment Effects

The heterogeneous treatment effect estimation procedure 
proposed by Imai and Ratkovic (2013)<DOI: 10.1214/12-AOAS593>.  
The proposed method is applicable, for example, when selecting 
a small number of most (or least) efficacious treatments from 
a large number of alternative treatments as well as when 
identifying subsets of the population who benefit (or are harmed 
by) a treatment of interest. The method adapts the Support Vector 
Machine classifier by placing separate LASSO constraints over the
pre-treatment parameters and causal heterogeneity parameters of
interest. This allows for the qualitative distinction between
causal and other parameters, thereby making the variable
selection suitable for the exploration of causal heterogeneity. The 
package also contains the function, CausalANOVA, which estimates 
the average marginal interaction effects by a regularized ANOVA as 
proposed by Egami and Imai (2016+). 
