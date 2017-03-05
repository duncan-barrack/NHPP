# Classification and clustering for samples of event time data using non-homogeneous Poisson process models
This repository contains Matab code accompaining the paper ... . If you use the code in your work please cite this paper. 

# Software and version numbers
The code has been tested on Matlab R2014a (see www.mathworks.com) but should work with subsequent versions. 

# Description of files
#### Classification ####
* `NHPP_train.m` - Obtains non homogeneous Poisson process (NHPP) rate function estimates for training set for classification tasks. It has the following dependencies:-
 * `NHPP_of.m` - objective function to be optmimised
 * `NHPP_con.m` - constraint for the optimisation
* `NHPP_test.m` - Obtains posteriori probabilities ofr test set data  

#### Clustering ####
* `NHPP_cluster.m` - Obtains NHPP rate function estimates and membership probabilities for the unsupervised learning task.  It has the following dependencies:-
 * `NHPP_of_EM.m` - Objective function for the NHPP mixture model
 * `NHPP_con_EM.m` - Optmisiation constraint
 
#### Synthetic data examples ####
* `classification_eg1.m`, `classification_eg2.m` and `classification_eg1.m` reproduce the classification results for synthetic data sets 1,2 and 3 resepctively from the paper. Similarly, `clustering_eg1.m`, `clustering_eg2.m` and `clustering_eg1.m` reproducde the clustering result for the synthetic data. These scripts have the following dependency:-
 * `NHPP.m` - Generates event times for a given rate function.

# Running the examples
To run `classification_eg1.m`, in Matlab, type 

`classification_eg1`

You should see the following

```
Fitting training set for Class 1
Running optimisation to obtain B-spline coefficients.
                                            First-order      Norm of
 Iter F-count            f(x)  Feasibility   optimality         step
    0     101    1.255458e+02    0.000e+00    1.348e+02
    1     202   -1.955777e+04    0.000e+00    8.037e-01    3.821e+02
    2     303   -1.957986e+04    0.000e+00    7.938e-01    4.048e+00
    3     404   -1.966988e+04    0.000e+00    7.434e-01    1.926e+01
    4     505   -1.979897e+04    0.000e+00    6.905e-01    4.254e+01
```
Once the optimisation is complete, the following raster plots showing a selection of event times for each sample as well as plots of the NHPP rate functions for each class will be produced.

![class_eg1_raster1_train](https://cloud.githubusercontent.com/assets/9549001/23590660/09436f0a-01dc-11e7-9d43-05b47a256164.jpg)


The other examples can be run in a similar way.
