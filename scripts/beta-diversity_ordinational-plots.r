library(tidyverse)
library(speedyseq)
# load data
load("data/fernow-analysis-data.RData")
source("scripts/ordinational_plot_functions.r")
#################################################
##### Analysis of bacterial beta diversity  #####
#################################################
permanova(ps, "site * sym * soil * date")
# site; sym; date are all significant factors.
# site:sym; site:date; sym:date are all significant interactions.
# site:sym:date is a marginally significant 3-way interaction.

permanova(ps %>% filter_sample_data(date == '2017'), "site * sym * soil")
# site and sym are both significant factors in bacterial community structure in 2017 samples
# additionally, site:sym is a significant interaction in 2017 samples
permanova(ps %>% filter_sample_data(date == '2022'), "site * sym * soil")
# site is a significant factor in bacterial community structure in 2022 samples
# no significant interactions are found in 2022 samples

# breaking down the samples by mutualist type type to evaluate whether site and soil type are significant factors in bacterial community structure
permanova(ps %>% filter_sample_data(date == '2017' & sym == 'AM' & soil == "Rhizosphere"), "site") %>% data.frame()
permanova(ps %>% filter_sample_data(date == '2017' & sym == 'AM' & soil == "Bulk"), "site") %>% data.frame()
permanova(ps %>% filter_sample_data(date == '2017' & sym == 'ECM'), "site * soil")
permanova(ps %>% filter_sample_data(date == '2022' & sym == 'AM'), "site * soil")
permanova(ps %>% filter_sample_data(date == '2022' & sym == 'ECM'), "site * soil")
# In all conditions tested above, site is a significant factor in bacterial community structure.

# breaking down the samples by sampling site to evaluate whether mutualist type and soil type are significant factors in bacterial community structure
permanova(ps %>% filter_sample_data(date == '2017' & site == 'ws7'), "sym * soil")
permanova(ps %>% filter_sample_data(date == '2017' & site != 'ws7'), "sym * soil")
permanova(ps %>% filter_sample_data(date == '2022' & site == 'ws7'), "sym * soil")
permanova(ps %>% filter_sample_data(date == '2022' & site != 'ws7'), "sym * soil")
# Mutualist type is a significant factor at both sites in 2017 samples, but is not significant in either sampling sites in 2022 samples.

#################################################################################################################################################################
#####   Ordinational plots of beta-diversity   ##################################################################################################################
#################################################################################################################################################################
plot.nmds(ps %>% filter_sample_data(date == '2017' & sym == 'AM'  & soil == "Rhizosphere"), "site", 'site')
plot.nmds(ps %>% filter_sample_data(date == '2017' & sym == 'AM'  & soil == "Bulk"), "site", 'site')

plot.nmds(ps %>% filter_sample_data(date == '2017' & sym == 'ECM'  & soil == "Rhizosphere"), "site", 'site')
plot.nmds(ps %>% filter_sample_data(date == '2017' & sym == 'ECM'  & soil == "Bulk"), "site", 'site')

plot.nmds(ps %>% filter_sample_data(date == '2022' & sym == 'AM'  & soil == "Rhizosphere"), "site", 'site')
plot.nmds(ps %>% filter_sample_data(date == '2022' & sym == 'AM'  & soil == "Bulk"), "site", 'site')

plot.nmds(ps %>% filter_sample_data(date == '2022' & sym == 'ECM'  & soil == "Rhizosphere"), "site", 'site')
plot.nmds(ps %>% filter_sample_data(date == '2022' & sym == 'ECM'  & soil == "Bulk"), "site", 'site')

site.1 = plot.nmds(ps %>% filter_sample_data(date == '2017' & sym == 'AM'  & soil == "Rhizosphere"), "site", 'site') + theme(legend.position = 'none')
site.2 = plot.nmds(ps %>% filter_sample_data(date == '2017' & sym == 'ECM'), "site * soil", 'site') + theme(legend.position = 'none')
site.3 = plot.nmds(ps %>% filter_sample_data(date == '2022' & sym == 'AM'), "site * soil", 'site') + theme(legend.position = 'none')
site.4 = plot.nmds(ps %>% filter_sample_data(date == '2022' & sym == 'ECM'), "site * soil", 'site') + theme(legend.position = 'none')
cowplot::plot_grid(site.1, site.2, site.3, site.4, nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"),
                        label_fontfamily = "Times", label_colour = "black", label_fontface = "plain")
ggsave("figures/ordinational_plots_site.svg", width = 6, height = 6, units = "in", dpi = 600)



plot.nmds(ps %>% filter_sample_data(date == '2017' & site == 'ws7'  & soil == "Rhizosphere"), "sym", 'sym')
plot.nmds(ps %>% filter_sample_data(date == '2017' & site == 'ws7'  & soil == "Bulk"), "sym", 'sym')

plot.nmds(ps %>% filter_sample_data(date == '2017' & site == 'ws3'  & soil == "Rhizosphere"), "sym", 'sym')
plot.nmds(ps %>% filter_sample_data(date == '2017' & site == 'ws3'  & soil == "Bulk"), "sym", 'sym')

plot.nmds(ps %>% filter_sample_data(date == '2022' & site == 'ws7'  & soil == "Rhizosphere"), "sym", 'sym')
plot.nmds(ps %>% filter_sample_data(date == '2022' & site == 'ws7'  & soil == "Bulk"), "sym", 'sym')

plot.nmds(ps %>% filter_sample_data(date == '2022' & site == 'ws3'  & soil == "Rhizosphere"), "sym", 'sym')
plot.nmds(ps %>% filter_sample_data(date == '2022' & site == 'ws3'  & soil == "Bulk"), "sym", 'sym')

sym.1 = plot.nmds(ps %>% filter_sample_data(date == '2017' & site == 'ws7'), "sym * soil", 'sym') + theme(legend.position = 'none')
sym.2 = plot.nmds(ps %>% filter_sample_data(date == '2017' & site == 'ws3'), "sym * soil", 'sym') + theme(legend.position = 'none')
sym.3 = plot.nmds(ps %>% filter_sample_data(date != '2017' & site == 'ws7'), "sym * soil", 'sym') + theme(legend.position = 'none')
sym.4 = plot.nmds(ps %>% filter_sample_data(date != '2017' & site == 'ws3'), "sym * soil", 'sym') + theme(legend.position = 'none')
cowplot::plot_grid(sym.1, sym.2, sym.3, sym.4, nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"),
                        label_fontfamily = "Times", label_colour = "black", label_fontface = "plain")
ggsave("figures/ordinational_plots_sym.svg", width = 6, height = 6, units = "in", dpi = 600)
