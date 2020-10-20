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
# The p values of Burden, KM, IEHC-Fisher, IEHC-adapt, IEHC-optim and ACAT

source("IEHC.R")

y <- read.table("data.txt",head=T)

X <- read.table("X.txt",head=T)

G <- read.table("G.txt",head=T)

M <- read.table("M.txt",head=T)

X = as.matrix(X)

G = as.matrix(G)

M = as.matrix(M)

fit = IEHC (data = y, X, G, M, combination_preference = "All")

$pvalue

  pvalue.Burden      3.202048e-06
  

  pvalue.KM           1.066000e-01
  

  pvalue.IEHC_optim   1.060116e-05
  

  pvalue.IEHC_adapt   6.461405e-06


  pvalue.IEHC_Fisher  1.003343e-05


  pvalue.ACAT         7.563492e-06
