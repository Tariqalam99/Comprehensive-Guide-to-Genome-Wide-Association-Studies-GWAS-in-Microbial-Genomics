# Comprehensive-Guide-to-Genome-Wide-Association-Studies-GWAS-in-Microbial-Genomics

Introduction

Welcome to this comprehensive guide on performing Genome-Wide Association Studies (GWAS) in microbial organisms. This tutorial will walk you through the installation of necessary packages, data loading, and running various statistical models to identify significant genetic associations. Whether you're a beginner or looking to refine your GWAS workflow, this guide provides step-by-step instructions to help you achieve accurate and meaningful results.

Prerequisites

Before you begin, ensure you have the following:

R installed on your computer. You can download it from CRAN.
Basic understanding of R programming and genomic data.
Access to your phenotypic and genotypic data files (e.g., FG95195.txt and Pnodorum_geno.txt).

1. Installing and Loading Required Packages

This demonstration primarily utilizes the GAPIT package for GWAS analysis. Additionally, we'll use the qqman package for creating Manhattan plots.

GAPIT Installation
GAPIT can be installed using two methods:

Option 1: Direct Sourcing (Easiest but requires sourcing every session)

source("http://www.zzlab.net/GAPIT/GAPIT.library.R")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")

Option 2: Using devtools (Permanent installation)

install.packages("devtools")
devtools::install_github("jiabowang/GAPIT3", force = TRUE)
library(GAPIT3)

Installing and Loading qqman
The qqman package is excellent for creating customizable Manhattan plots.

install.packages("qqman")
library(qqman)

2. Loading and Viewing Data

Ensure your working directory is set to the location of your data files.

# Set working directory (modify the path as needed)
setwd("/path/to/your/data")

# Load phenotypic and genotypic data
myY <- read.table("FG95195.txt", header = TRUE)
myG <- read.table("Pnodorum_geno.txt", header = FALSE)

# View the datasets
View(myY)
View(myG)

3. Running a General Linear Model (GLM)

GWAS typically involves testing multiple statistical models. We'll start with a naive GLM that does not account for population structure or kinship.

glm_naive <- GAPIT(Y = myY, G = myG, PCA.total = 0, model = "GLM")

Analyzing Results
Focus on QQ plots and Manhattan plots to assess associations.

# QQ Plot
qq(glm_naive$GWAS$P.value)

# Convert chromosome and position to numeric
glm_naive$GWAS$Chromosome <- as.numeric(glm_naive$GWAS$Chromosome)
glm_naive$GWAS$Position <- as.numeric(glm_naive$GWAS$Position)

# Manhattan Plot
manhattan(glm_naive$GWAS, chr = "Chromosome", bp = "Position", snp = "SNP", p = "P.value", genomewideline = 6.8)

Interpretation: The QQ plot compares observed p-values to expected ones under the null hypothesis. Deviations indicate potential associations. The Manhattan plot visualizes significant markers across chromosomes.

4. Incorporating Population Structure with GLM

To correct for potential population stratification, include Principal Components (PCs) as covariates.

glm_model_selection <- GAPIT(Y = myY, G = myG, PCA.total = 3, Model.selection = TRUE)

Comparing Models

# QQ Plot for Model Selection
qq(glm_model_selection$GWAS$P.value)

# Convert chromosome and position to numeric
glm_model_selection$GWAS$Chromosome <- as.numeric(glm_model_selection$GWAS$Chromosome)
glm_model_selection$GWAS$Position <- as.numeric(glm_model_selection$GWAS$Position)

# Manhattan Plot for Model Selection
manhattan(glm_model_selection$GWAS, chr = "Chromosome", bp = "Position", snp = "SNP", p = "P.value", genomewideline = 6.8)

Interpretation: Including PCs helps control for false positives by accounting for population structure, leading to more reliable associations.

5. Running a Mixed Linear Model (MLM)

To further control for relatedness (kinship) among samples, use a Mixed Linear Model.

mlm_pca <- GAPIT(Y = myY, G = myG, PCA.total = 1, model = "MLM", kinship.algorithm = "EMMA")

Comparing All Models

# QQ Plot for MLM
qq(mlm_pca$GWAS$P.value)

# Convert chromosome and position to numeric
mlm_pca$GWAS$Chromosome <- as.numeric(glm_naive$GWAS$Chromosome)
mlm_pca$GWAS$Position <- as.numeric(glm_naive$GWAS$Position)

# Manhattan Plot for MLM
manhattan(mlm_pca$GWAS, chr = "Chromosome", bp = "Position", snp = "SNP", p = "P.value", genomewideline = 6.8)

Interpretation: The MLM often provides the best control for false positives by accounting for both population structure and kinship, leading to more accurate identification of significant loci.

Conclusion

Through this guide, you've learned how to perform GWAS in microbial organisms using different statistical models to identify significant genetic associations. By comparing GLM and MLM approaches, you can select the most appropriate model for your specific dataset and research objectives.
