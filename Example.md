# Example of IEHC
# Input data
IEHC requires data (survival time and status, given in data.txt), X (covariates, X should be a n by p matrix), G (the weighted burden score, G should be a n by 1 matrix) and M (genotypes). For the input data, no missing data is allowed. So, missing data should be removed before data analysis.
# The function of IEHC
IEHC = function (

  data, X, G, M,
  
  chisq_app = "3M",
  
  combinations_return=TRUE,
  
  combination_preference = "All", # other options: "optim", "adapt", "Fisher", "Burden", "KM", "ACAT"
  
  acc = 5e-10,
  
  accurate_app_threshold = -log10(0.05),
  
  max_core = 4
){
