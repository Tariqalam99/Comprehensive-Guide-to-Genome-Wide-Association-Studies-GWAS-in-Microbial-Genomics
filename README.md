# Comprehensive-Guide-to-Genome-Wide-Association-Studies-GWAS-in-Microbial-Genomics

# Introduction

Welcome to this comprehensive guide on performing Genome-Wide Association Studies (GWAS) in microbial organisms. This tutorial will walk you through the installation of necessary packages, data loading, and running various statistical models to identify significant genetic associations. Whether you're a beginner or looking to refine your GWAS workflow, this guide provides step-by-step instructions to help you achieve accurate and meaningful results.

## Prerequisites

Before you begin, ensure you have the following:

- **R** installed on your computer. You can download it from [CRAN](https://cran.r-project.org/).
- Basic understanding of R programming and genomic data.
- Access to your phenotypic and genotypic data files (e.g., `FG95195.txt` and `Pnodorum_geno.txt`).


## Installing/loading packages

This demonstration will primarily use the package GAPIT. The user manual is very descriptive and easy to follow (https://www.zzlab.net/GAPIT/gapit_help_document.pdf)


Lipka, AE., Tian, F., Wang, Q., Peiffer, J., Li, M., Bradbury, P.J., Gore, M.A., Buckler, E.S., and Zhang, Z. GAPIT: genome association and prediction integrated tool. Bioinformatics. 2012. 28(18):2397-2399.

Install/load GAPIT with one of the following options (first option is easiest, but you have to do it every time you start a new R session):

Option 1:

```{r, eval = FALSE}
source("http://www.zzlab.net/GAPIT/GAPIT.library.R")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
```

Option 2:

```{r, eval = FALSE}
install.packages("devtools")
devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(GAPIT3)
```

```{r, eval=TRUE, echo=FALSE, message = FALSE, error=FALSE, warning=FALSE}
library(GAPIT3)
```

The 'qqman' package has a nice function for making Manhattan plots. GAPIT will automatically output Manhattan plots, but I find it a little bit easier to produce customizable figures using this package.

```{r, eval = FALSE}
install.packages("qqman")
library(qqman)
```

```{r, echo = FALSE, message = FALSE, error=FALSE, warning=FALSE}
library(qqman)
```
## Loading and viewing data

First, let's load our phenotypic and genotypic data (remember to set your working directory to where the data is located)

```{r, eval = TRUE, cache = TRUE}
myY <- read.table("FG95195.txt", head = TRUE)
myG <- read.table("Pnodorum_geno.txt", head = FALSE)
```


Next, we will examine each dataset to get familiar with the format

```{r, eval = FALSE}
View(myY)
View(myG)
```

## Running a general linear model (GLM)

A typical GWAS workflow includes testing various models. First, we will test two different general linear models. One will only test associations between the markers and traits, without correction for population structure or kinship. This is commonly referred to as a naive model.

```{r, eval = TRUE, results = "hide", cache = TRUE, error = FALSE, message = FALSE, warning = FALSE}
glm_naive <- GAPIT(Y=myY, G=myG, PCA.total = 0, model = "GLM")
```

This analysis will output a lot of useful files and figures. For the purposes of this workshop, we will focus on the QQ plot and Manhattan plots.

```{r, echo=TRUE, cache = TRUE}
qq(glm_naive$GWAS$P.value)
```

This QQ plot is comparing the observed p-values from this analysis to the expected distribution of p-values under the null hypothesis. The solid line represents a p-value distribution if no significant associations were detected. A deviation above that line would represent over-inflation of p-values, indicating the potential detection of false positives. If we identify a true association, we do expect to see a sharp peak towards in the upper right of this figure. Overall, based on this plot, we see that we do have over-inflated p-values and we should try to correct this using a different model.

```{r, echo=TRUE, cache = TRUE}
glm_naive$GWAS$Chromosome <- as.numeric(glm_naive$GWAS$Chromosome)
glm_naive$GWAS$Position <- as.numeric(glm_naive$GWAS$Position)
```

```{r, cache=TRUE, message=FALSE,error=FALSE,warning=FALSE}
manhattan(glm_naive$GWAS, chr="Chromosome", bp="Position", snp="SNP", p = "P.value", genomewideline = 6.8)
```

The Manhattan plot shows us our significant associations. The chromosomes/contigs are listed on the x-axis and the -log10(p) values are listed on the y-axis. Each dot represents a single marker. We can see that we identified a significant association in the sub-telomeric region of chromosome 3, as well as some other loci that crossed our significance threshold (Bonferroni correction).


Since we could see that our p-values may be inflated, the next GLM will incorporate population structure as a covariate (fixed effect). Principal coordinates analysis (PCA) is a common method to do this and is easily incorporated in GAPIT. Other options could be output from programs such as STRUCTURE, which could also be imported into R and used as a covariate.

Since we can use results from PCA as a covariate, how many PCs should we decide to use? Each PC explains a proportion of the cumulative variation. Typically, you could select the amount of PCs that explain ~25% or 50% of the variation. In GAPIT, we can activate 'model selection', which will select the 'best' model based on Bayesian information criterion (BIC).

```{r, eval = TRUE, results = "hide", cache = TRUE, error = FALSE, message = FALSE, warning = FALSE}
glm_model_selection <- GAPIT(Y=myY, G=myG, PCA.total = 3, Model.selection = TRUE)
```

Now, compare our results from the naive model to the GLM with correction for population structure.

```{r, echo=TRUE, cache = TRUE}
qq(glm_model_selection$GWAS$P.value)
```

It looks like this model did a better job at controlling inflated p-values (i.e. potential spurious associations). However, there is still some deviation from the null distribution.


```{r, echo=TRUE, cache = TRUE}
glm_model_selection$GWAS$Chromosome <- as.numeric(glm_naive$GWAS$Chromosome)
glm_model_selection$GWAS$Position <- as.numeric(glm_naive$GWAS$Position)
```

```{r, cache=TRUE, message=FALSE,error=FALSE,warning=FALSE}
manhattan(glm_model_selection$GWAS, chr="Chromosome", bp="Position", snp="SNP", p = "P.value", genomewideline = 6.8)
```

Our major peak on chromosome 3 is still highly significant. Some of the minor association detected with the naive model are now below the significance threshold and may have been spurious associations.

## Running a mixed linear model (MLM)

In addition to incorporating population structure as a fixed effect, we can also control for false positives by including kinship (relatedness) as a random effect in a mixed linear model. There are several variations of this MLM (compressed MLM, enhanced MLM) that can be run in GAPIT, but we will just focus on the 'classic' example using the EMMA algorithm.

```{r, eval = TRUE, results = "hide", cache = TRUE, error = FALSE, message = FALSE, warning = FALSE}
mlm_pca <- GAPIT(Y=myY, G=myG, PCA.total = 1, model = "MLM", kinship.algorithm = "EMMA")
```

Compare these results to both iterations of the GLM. What differences do you observe in the QQ plots? What differences do you observe in the Manhattan plots? Is there a 'best' model?

```{r, echo=TRUE, cache = TRUE}
qq(mlm_pca$GWAS$P.value)
```

```{r, echo=TRUE, cache = TRUE}
mlm_pca$GWAS$Chromosome <- as.numeric(glm_naive$GWAS$Chromosome)
mlm_pca$GWAS$Position <- as.numeric(glm_naive$GWAS$Position)
```

```{r, cache=TRUE, message=FALSE,error=FALSE,warning=FALSE}
manhattan(mlm_pca$GWAS, chr="Chromosome", bp="Position", snp="SNP", p = "P.value", genomewideline = 6.8)
```

After testing all three models, it appears that the MLM may be the 'best' model. But what does that really mean? All three models detected our major locus. From a biological perspective, if we are aiming to identify the major effect genes underlying virulence, any of these models would have resulted in the successful detection of this gene. However, if we wanted to investigate the minor effect loci further, it may be best to start with the model that appeared to control for false positives more effectively (but still consider all models!).
