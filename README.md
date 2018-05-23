# Sparse Estimation for Functional Semiparametric Additive Models

This repository consists of coding for two simulation studies and two applications in the manuscript.

1. sim1.R: the R code for the simulation study for functional semiparametric additive model with scalar predictors;

2. sim2.R: the R code for the simulation study for functional semiparametric additive model without scalar predictors;

3. tecator.txt: the tecator data. The description for the tecator data is available on http://lib.stat.cmu.edu/datasets/tecator;

4. The folder "MatlabCode" contains one matlab file (tecator_script_tune.m), which is used to implment PACE for both training data and test data. This folder also contains a folder "pace_release2.11", which is the PACE package.

5. trainRes.mat: results from matlab codes for the training data.

6. testpace.mat: results from matlab codes for the test data.

7. Tecator analysis.R: the R code for analyzing the tecator data. We also compared different estimation methods.

8. ADHD analysis.R: the R code for analyzing the ADHD data. We also compared different estimation methods.

9. Gaussian.R, internal gau.r and iterative update.R: These three files are some R functions used for our proposed method, which is called FSAM-COSSO in our manuscript. 
