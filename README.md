# Classification and clustering for samples of event time data using non-homogeneous Poisson process models
This repository contains Matab code for the paper "Classification and clustering for samples of event time data using non-homogeneous Poisson process models" http://arxiv.org/abs/1703.02111. If you use the code in your work please cite this paper. 

# Software and version numbers
The code has been tested on Matlab R2014a (see www.mathworks.com) but should work with subsequent versions. 

# Description of files
#### Classification ####
* `NHPP_train.m` - Obtains non homogeneous Poisson process (NHPP) rate function estimates for training set data. It has the following dependencies:-
 * `NHPP_of.m` - objective function to be optimised.
 * `NHPP_con.m` - constraint for the optimisation.
* `NHPP_test.m` - Obtains posteriori probabilities for test set data.  

#### Clustering ####
* `NHPP_cluster.m` - Obtains NHPP rate function estimates and membership probabilities for each sample for the unsupervised learning task.  It has the following dependencies:-
 * `NHPP_of_EM.m` - Objective function for the NHPP mixture model.
 * `NHPP_con_EM.m` - Optimisation constraint.
 
#### Synthetic data examples ####
* `classification_eg1.m`, `classification_eg2.m` and `classification_eg1.m` reproduce the classification results for synthetic data sets 1,2 and 3 respectively from the paper. Similarly, `clustering_eg1.m`, `clustering_eg2.m` and `clustering_eg1.m` reproduce the clustering results for the synthetic data. These scripts have the following dependency:-
 * `NHPP.m` - Generates event times for a given rate function.

# Running the examples
To run `classification_eg1.m`, in Matlab, type 

`classification_eg1`

You should see the following.

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

![class_eg1_raster1_train](https://cloud.githubusercontent.com/assets/9549001/23590723/392b8fda-01dd-11e7-87dd-cf68c3bcff9c.jpg) ![class_eg1_raster2_train](https://cloud.githubusercontent.com/assets/9549001/23590932/9ff992f4-01e0-11e7-8fe3-a4911869897c.jpg)

![class_eg1_raster1_test](https://cloud.githubusercontent.com/assets/9549001/23590935/a46d9f24-01e0-11e7-87b0-0a01d999e19b.jpg) ![class_eg1_raster2_test](https://cloud.githubusercontent.com/assets/9549001/23590936/a9d90ae8-01e0-11e7-9916-44676f0a8b3d.jpg)

![class_eg1_nhpp1](https://cloud.githubusercontent.com/assets/9549001/23590939/aeb4333a-01e0-11e7-9d32-b61778220ed6.jpg) ![class_eg1_nhpp2](https://cloud.githubusercontent.com/assets/9549001/23590942/b41dca52-01e0-11e7-9fca-a20281c4470d.jpg)

The other examples can be run in a similar way.
