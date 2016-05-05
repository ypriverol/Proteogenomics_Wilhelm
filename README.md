# Proteogenomics_Wilhelm

-------

# FINAL FOLDER:

The analysis used in the paper. 

DATA PREPARATION

1_ Get the data from their files and put it in the format they used in their analysis

2_ Normalize rows of the data to remove the differences in expression between proteins

PREDICTION

3_ Prediction methods (ratio, ratio loo, control, control loo)

EVALUATION

4_ Evaluations ACROSS genes (Correlation)

5_ Example plot from real data

6_ Example plots from simulated data

7_ Full model (one model across genes & tissues) 

8_ Evaluation PER gene (correlation, R2, MSE,...), also between vs per-gene variance

9_ Analysis of the variation of ratios (per gene across tissues)

10_ Predicting protein from random genes - to predict protein of gene X, averaging across all mRNAs from genes other than X

11_ Predicting protein from random genes - using just one mRNA from a gene other than X to predict protein of gene X

12_ Using the randomized predictions above to derive a p-value for the correlations we do find


# OLD FOLDER:

Here I was trying different methods to predict protein from mRNA quantity. This part is more experimental

Focus on files starting with 4. They summarize results.

The file 0_Pipeline.R creates a knitr report from all other files in this repository.

Other files that start with a "0" are old files that try to reproduce the origin results in the publication and explore effects of variability and missing values on the result.

-------

1_ are preping data for prediction; here only the genes are used where there are no missing values in mRNA or protein!

2_ are making predictions

3_ are evaluating predictions

4_ Summary of the results. Here, all genes are used, as it was originally done by Wilhelm et al. Effects of missing values and variability are also examined.
