---
title: "Fernow Forest Microbiome Data Analysis"
author: "Michael D. Millican"
output: beamer_presentation
bibliography: citations/fernow-citations.bib
---

```{r library_load, include=FALSE}
library(tidyverse)
library(vegan)
library(speedyseq)
```
# Required R packages
- The following packages are required to run the following scripts
    - devtools
        - needed to install speedyseq package
    - BiocManager
        - needed to install phyloseq package
    - tidyverse
    - vegan
    - phyloseq
    - speedyseq
## Installing packages
- The following script will install R packages
```
{r install_packages, eval=FALSE}
    install.packages(c("devtools", "BiocManager"))
    install.packages("tidyverse")
    install.packages("vegan")
    BiocManager::install("phyloseq")
    devtools::install_github("mikemc/speedyseq")
```
---
# Data categories
Key for categorical data codes:
site = which watershed did sample come from
sym = dominate mycorrhizal fungi expected in sampling location
soil = soil bioshpere sampled (rhizosphere or bulk soil)
date = sampling year (?)
```{r loading_data, echo=FLASE, include=FALSE}
load("data/fernow-analysis-data.RData")
```
---
# Alpha Diversity
### Statistical analysis and data visulaization of bacterial alpha diversity

---
# Alpha Diversity Stats
---
### Bacterial richness
Does bacterial richness vary significantly based on the combined effects of sampling site, symbiotroph, soil biosphere, and sampling date, and if so, which factors or interactions are driving these variations?

```{r alpha_stats, echo=FALSE,warning=FALSE,error=FALSE}
df.meta %>% 
    lm(rich ~ site * sym * soil * date, data = .) %>%
    car::Anova()
```
- sym, soil are both significant for bacterial richness (p < 0.05)
- site:date and sym:date are both significant interactions for richness.
---
### Is there a significant difference in bacterial richness among the different sampling sites for 2017 sampling date?
```{r alpha_anova_2017, echo=FALSE, warning=FALSE}
df.meta %>% filter(date == '2017') %>%
    lm(rich ~ site * sym * soil, data = .) %>%
    car::Anova()
```
---
### Is there a significant difference in bacterial richness among the different sampling sites for 2017 sampling date?
```{r alpha_anova_2017, echo=FALSE, warning=FALSE}
df.meta %>% filter(date != '2017') %>%
    lm(rich ~ site * sym * soil, data = .) %>%
    car::Anova()
```
---
### Bacterial Diversity
For bacterial alpha diversity we are using the estimated number of species from the probability of interspecific encounter, or ENSpie, metric (also known as inverse simpson's).

For more reading on ENSpie see the following references:
[Chase, J and Knight, T. 2013. Scale-dependent effect sizes of ecological drivers on biodiversity: why standardised sampling is not enough](https://doi.org/10.1111/ele.12112)
[Jost, L. 2006. Entropy and diversity](https://doi.org/10.1111/j.2006.0030-1299.14714.x)
[Seabloom *et al.* 2023. Globally consistent response of plant microbiome diversity across hosts and continents to soil nutrients and herbivores.](https://doi.org/10.1038/s41467-023-39179-w)
---
Does bacterial diversity vary significantly based on the combined effects of sampling site, symbiotroph, soil biosphere, and sampling date, and if so, which factors or interactions are driving these variations?
```{r alpha_stats, echo=FALSE,warning=FALSE,error=FALSE}
df.meta %>% 
    lm(ens ~ site * sym * soil * date, data = .) %>%
    car::Anova()
```
- No categorical variables are found to be significant for bacterial diversity.
- site, sym:soil, site:date were all marginally significant (p < 0.1; p > 0.05)
---
### Is there a significant difference in bacterial diversity among the different sampling sites for 2017 sampling date?
```{r alpha_anova_2017, echo=FALSE, warning=FALSE}
df.meta %>% filter(date == '2017') %>%
    lm(ens ~ site * sym * soil, data = .) %>%
    car::Anova()
```
---
### Is there a significant difference in bacterial diversity among the different sampling sites for 2017 sampling date?
```{r alpha_anova_2017, echo=FALSE, warning=FALSE}
df.meta %>% filter(date != '2017') %>%
    lm(ens ~ site * sym * soil, data = .) %>%
    car::Anova()
```
---
### Breaking the data down by sampling date and soil type:
```{r, echo=FALSE, warning=FALSE}
df.meta %>% filter(date == '2017' & soil == 'Bulk') %>%
    lm(ens ~ site * sym, data = .) %>%
    car::Anova()
```
---
```{r, echo=FALSE, warning=FALSE}
df.meta %>% filter(date == '2017' & soil == 'Rhizosphere') %>%
    lm(ens ~ site * sym, data = .) %>%
    car::Anova()
```
---
```{r, echo=FALSE, warning=FALSE}
df.meta %>% filter(date != '2017' & soil == 'Bulk') %>%
    lm(ens ~ site * sym, data = .) %>%
    car::Anova()
```
---
```{r, echo=FALSE, warning=FALSE}
df.meta %>% filter(date != '2017' & soil == 'Rhizosphere') %>%
    lm(ens ~ site * sym, data = .) %>%
    car::Anova()
```
---
### Sampling site is only significant for bulk soil in 2022 samples.

