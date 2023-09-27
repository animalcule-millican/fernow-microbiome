# load packages
library(tidyverse)
library(speedyseq)
# load data
load("data/fernow-analysis-data.RData")
#################################################
##### Analysis of bacterial Alpha diversity #####
#################################################
# Influence of experimental factors on bacterial diversity
df.meta %>% 
    lm(rich ~ site * sym * soil * date, data = .) %>%
    car::Anova()
## sym, soil are both significant for bacterial richness (p < 0.05)
## site:date and sym:date are both significant interactions for richness.
## No marginally significant categorical variables are found for richness.
df.meta %>% filter(date == '2017') %>%
    lm(rich ~ site * sym * soil, data = .) %>%
    car::Anova()
## sym, soil are both significant for bacterial richness (p < 0.05)
df.meta %>% filter(date != '2017') %>%
    lm(rich ~ site * sym * soil, data = .) %>%
    car::Anova()
## soil is marginally significant for bacterial richness (p < 0.1; p > 0.05)
#########
## Conclusions:
## For the 2017 samples, the type of mutualist type the community was associated with and the watershed site both have a significant effect on bacterial richness. These two categorical factors do not interact suggesting independant influence on the bacterial richness.
## For the 2022 samples, whether samples came from a rhizosphere site or a bulk soil site has a marginal effect on bacterial richness.

# Influence of experimental factors on bacterial diversity
df.meta %>% 
    lm(ens ~ site * sym * date * soil, data = .) %>%
    car::Anova()
## No categorical variables are found to be significant for bacterial diversity.
## site, sym:soil, site:date were all marginally significant (p < 0.1; p > 0.05)
df.meta %>% filter(date == '2017' & soil == 'Bulk') %>%
    lm(ens ~ site * sym, data = .) %>%
    car::Anova()
df.meta %>% filter(date == '2017' & soil == 'Rhizosphere') %>%
    lm(ens ~ site * sym, data = .) %>%
    car::Anova()
df.meta %>% filter(date != '2017' & soil == 'Bulk') %>%
    lm(ens ~ site * sym, data = .) %>%
    car::Anova()
df.meta %>% filter(date != '2017' & soil == 'Rhizosphere') %>%
    lm(ens ~ site * sym, data = .) %>%
    car::Anova()
##############################################################################################################################
##############################################################################################################################
#####   Diversity Figures   ##################################################################################################
##############################################################################################################################
##############################################################################################################################
df.meta %>% 
    group_by(site, sym, date, soil) %>%
    summarise(mean = mean(rich), se = sd(rich)/sqrt(n())) %>%
    ggplot(aes(x = site, y = mean, fill = sym)) +
    geom_linerange(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(width = 0.5), linewidth = 1) +
    geom_point(position = position_dodge(width = 0.5), size = 5, shape = 21) +
    facet_grid(date~soil) + ylab("Observed richness") + theme(axis.title.x = element_blank()) + theme_bw()
ggsave("figures/observed-richness.svg", width = 3, height = 3, units = "in", dpi = 600)

df.meta %>% 
    group_by(site, sym, date, soil) %>%
    summarise(mean = mean(ens), se = sd(ens)/sqrt(n())) %>%
    ggplot(aes(x = site, y = mean, fill = sym)) +
    geom_linerange(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(width = 0.5), linewidth = 1) +
    geom_point(position = position_dodge(width = 0.5), size = 5, shape = 21) +
    facet_grid(date~soil) + ylab("Estimated number of species (ENSpie)") + theme(axis.title.x = element_blank()) + theme_bw()
ggsave("figures/ens-diversity.svg", width = 3, height = 3, units = "in", dpi = 600)

