# Initial differences in bacterial richness and diversity
library(tidyverse)
library(speedyseq)
# load data
load("data/fernow-analysis-data.RData")
"""
What are the initial differences in bacterial richness?
    - between soil types overall
    - between watersheds overall
    - between symbiotroph ecosystems overall
    - between soil types in each watershed
    - between soil types in each watershed for each symbiotroph ecosystem
    - between symbiotroph ecosystems in each watershed
"""
df.meta %>% filter(date == '2017') %>%
    lm(rich ~ soil * site * sym, data = .) %>%
    car::Anova()
"""
There is a significant difference in bacterial richness between watershed site and symbiotroph ecosystem, but not between soil types overall.
There are no significant interaction terms. 
"""

df.meta %>% filter(date == '2017') %>%
    lm(rich ~ site * sym * soil, data = .) %>%
    emmeans::emmeans(., ~ sym | site * soil) %>%
    multcomp::cld(., Letters = letters)
"""
There are significant differences in richness between symbiotroph ecosystems in the Bulk soil of watershed 7, but not watershed 3. 
Additionally, there are significant differences in richness between symbiotroph ecosystems in the Rhizosphere soil of watershed 3 and 7.
In all situations bacterial richness in ECM dominated ecosystems is higher than AMF dominated ecosystems.
"""
df.meta %>% filter(date == '2017') %>%
    lm(rich ~ site * sym * soil, data = .) %>%
    emmeans::emmeans(., ~ site | sym * soil) %>%
    multcomp::cld(., Letters = letters)
"""
A significant difference in richness betwee ws3 and ws7 for Bulk soil in AMF dominated ecosystems.
No other richness differences between sites were observed to be significant.
"""

########################################################################################################################
########################################################################################################################

df.meta %>% filter(date == '2017') %>%
    lm(ens ~ soil * site * sym, data = .) %>%
    car::Anova()
"""
There are no significant differences in bacterial diversity between soil types, site, symbiotroph ecosystem, or interactions between terms.
"""
