# Later differences in bacterial richness and diversity
library(tidyverse)
library(speedyseq)
# load data
load("data/fernow-analysis-data.RData")
"""
What are the later differences in bacterial richness?
    - between soil types overall
    - between watersheds overall
    - between symbiotroph ecosystems overall
    - between soil types in each watershed
    - between soil types in each watershed for each symbiotroph ecosystem
    - between symbiotroph ecosystems in each watershed
"""
df.meta %>% filter(date == '2022') %>%
    lm(rich ~ soil * site * sym, data = .) %>%
    car::Anova()

df.meta %>% filter(date == '2022') %>%
    lm(ens ~ soil * site * sym, data = .) %>%
    car::Anova()
"""
There is a marginal difference in bacterial richness between soil types overall.
There is a significant difference in bacterial ENSpie diversity between watersheds overall.
"""

df.meta %>% filter(date == '2022') %>%
    lm(rich ~ site * sym * soil, data = .) %>%
    emmeans::emmeans(., ~ soil | site * sym) %>%
    multcomp::cld(., Letters = letters)
"""
There are no significant pairwise differences between soil types at each site and symbiotroph factor level.
"""
df.meta %>% filter(date == '2022') %>%
    lm(ens ~ site * sym * soil, data = .) %>%
    emmeans::emmeans(., ~ site | soil * sym) %>%
    multcomp::cld(., Letters = letters)
"""
There are no significant pairwise differences between watersheds at each soil and symbiotroph factor level.
"""


"""
At the alpha diversity level, there are minimal differences in bacterial communities across the watersheds, soil types, and symbiotroph ecosystems.
"""