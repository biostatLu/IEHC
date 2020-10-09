# IEHC: integrative eQTL hierarchical Cox model for SNP set association based on mixed models
# Introduction
Integrative eQTL hierarchical Cox (IEHC) is a R procedure for examining the association of a set of genetic variants based on genes under the framework of mixed models. Very similar to the popular MiST method which conducts the burden and the variance components test to binary or continuous outcome. IEHC uses Integrative eQTL hierarchical Cox model to test for the association by testing both the burden and the variance components according to the survival risk outcome. 

Specifically, let y be a n by 1 vector of survival outcome on n individuals, X is a n by p matrix for covariates, G is a n by S matrix for genotypes of SNPs for a genetic region (i.e. gene), θ quantifies the association between the survival risk and the weighted burden score Gβ which explained by eQTL information, b quantifies the association between the survival risk and genotypes G which not explained by eQTL information and c a p-vector of fixed effect sizes for clinical covariates. We relate y, Gβ, X and G by a Cox mixed model:

y = (Gβ) × θ + Gb + Xc, b ~ N(0, τ)

Above, τ is the genetic variance which not explained eQTL information.

IEHC examines the association of G and Gβ with y (while controlling for X) by testing for:

H0: θ = 0 and b = 0 <==> H0: θ = 0 and τ = 0

This is a joint test including both fixed effect and random effects: the first component of H0 examines the influence of genetic variants on the survival risk explained by eQTLs; while the second component examines the impact of genetic variants beyond the effects of eQTLs. Briefly, we derive the test statistic for θ under H0: θ = 0 and τ = 0 as usual, while we derive the score statistic for τ under τ = 0 but without the constraint of θ = 0. By doing this, we ensure that these two statistics are independent. This strategy substantially eases the development of test statistics for the joint test. In conclusion, under this framework two asymptotically independent statistics can be derived: one for a scale (i.e. θ) in the general Cox model and the other for the variance component (i.e. τ) in the KM Cox model.

In order to aggregate the two independent test statistics, we propose three p-value combination approaches (i.e. IEHC-Fisher, IEHC-adapt and IEHC-optim).

# References
DR Cox. Regression Models and Life-Tables, Journal of the royal statistical society. Series B (Methodological), 1972, 34(2): 187-220. 

Tianxi Cai, Giulia Tonini and Xihong Lin. Kernel machine approach to testing the significance of multiple genetic markers for risk prediction. Biometrics, 2011, 67(3): 975-986.

Xinyi Lin, Tianxi Cai, Michael C. Wu, Qian Zhou, Geoffrey Liu, David C. Christiani and Xihong Lin. Kernel machine SNP-set analysis for censored survival outcomes in genome-wide association studies. Genetic Epidemiology, 2011, 35(7), 620-631.

Ping Zeng, Yang Zhao, Jin Liu, Liya Liu, Liwei Zhang, Ting Wang, Shuiping Huang and Feng Chen. Likelihood Ratio Tests in Rare Variant Detection for Continuous Phenotypes. Annals of Human Genetics, 2014, 78(5): 320-332. DOI: 10.1111/ahg.12071

Jianping Sun, Yingye Zheng, Li Hsu. A unified mixed-effects model for rare-variant association in sequencing studies. Genetic Epidemiology, 2013, 37(4), 334-344.

Yu-Ru Su, Chongzhi Di, Stephanie Bien, Licai Huang, Xinyuan Dong, Goncalo Abecasis, Sonja Berndt, Stephane Bezieau, Hermann Brenner, Bette Caan, Graham Casey, Jenny Chang-Claude, Stephen Chanock, Sai Chen, Charles Connolly, Keith Curtis, Jane Figueiredo, Manish Gala, Steven Gallinger, Tabitha Harrison, Michael Hoffmeister, John Hopper, Jeroen R. Huyghe, Mark Jenkins, Amit Joshi, Loic Le Marchand, Polly Newcomb, Deborah Nickerson, John Potter, Robert Schoen, Martha Slattery, Emily White, Brent Zanke, Ulrike Peters and Li Hsu. A Mixed-Effects Model for Powerful Association Tests in Integrative Functional Genomics. The American Journal of Human Genetics, 2018, 102(5), 904-919.
