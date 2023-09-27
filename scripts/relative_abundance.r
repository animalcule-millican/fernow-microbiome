library("tidyverse")
library("speedyseq")
load("data/fernow-analysis-data.RData")

###########################################
######### Bacterial Phyla #################
###########################################
ps %>% 
    speedyseq::filter_tax_table(!is.na(phylum)) %>%
    speedyseq::tax_glom(taxrank = "phylum") %>%
    speedyseq::transform_sample_counts(function(x) x / sum(x)) %>%
    speedyseq::psmelt() %>%
    group_by(site, sym, soil, date, phylum, Sample)%>%
    dplyr::summarize(abun=sum(Abundance))%>%
    group_by(site, sym, soil, date, phylum)%>%
    dplyr::summarize(mean.abun = mean(abun)) %>%
    data.frame %>% arrange(-mean.abun)
    #write.csv("data/relative_abundance_fungal_phylum.csv")
########################################
df.ab = ps %>% 
    speedyseq::filter_tax_table(!is.na(phylum)) %>%
    speedyseq::tax_glom(taxrank = "phylum") %>%
    speedyseq::transform_sample_counts(function(x) x / sum(x)) %>%
    speedyseq::psmelt() %>% filter(date == '2017' & soil == 'Rhizosphere') %>%
    data.frame %>% 
    group_by(site, sym, phylum) %>% 
    mutate(median =  median(Abundance)) %>%
    mutate(phylum = ifelse(median <= 0.01, 'Remainder', phylum)) %>%
    data.frame %>%
    group_by(site, sym, phylum, Sample)%>%
    dplyr::summarize(abun=sum(Abundance))%>%
    group_by(site, sym, phylum)%>%
    dplyr::summarize(mean.abun = mean(abun)) %>%
    data.frame %>% arrange(-mean.abun)
df.ab$phylum %>% unique()
phy = c(df.ab %>% filter(phylum != 'Remainder') %>% arrange(-mean.abun) %>% .$phylum %>% unique(), 'Remainder')

options(
  ggplot2.discrete.colour = pals::polychrome(25) %>% data.frame %>% .[3:25,1],
  ggplot2.discrete.fill = c("#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD", "#AA0DFE",'#353636')
)

df.ab %>% mutate(phylum = factor(phylum, levels = phy)) %>%
 ggplot(aes(x = site, y = mean.abun, fill = phylum)) +
        geom_bar(stat = 'identity', color = '#353636', alpha = 0.8) + 
        facet_grid(.~sym) + 
        ylab("Relative abundance") +
        guides(fill=guide_legend(title="phylum")) + 
        scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.1))) + theme_bw()
ggsave(file = 'figures/relative-abundance_bacteria_site-sym_2017-Rhizosphere_phylum.svg', width = 6, height = 3, dpi = 900, units = 'in')

########################################
df.ab = ps %>% 
    speedyseq::filter_tax_table(!is.na(phylum)) %>%
    speedyseq::tax_glom(taxrank = "phylum") %>%
    speedyseq::transform_sample_counts(function(x) x / sum(x)) %>%
    speedyseq::psmelt() %>% filter(date == '2017' & soil == 'Bulk') %>%
    data.frame %>% 
    group_by(site, sym, phylum) %>% 
    mutate(median =  median(Abundance)) %>%
    mutate(phylum = ifelse(median <= 0.01, 'Remainder', phylum)) %>%
    data.frame %>%
    group_by(site, sym, phylum, Sample)%>%
    dplyr::summarize(abun=sum(Abundance))%>%
    group_by(site, sym, phylum)%>%
    dplyr::summarize(mean.abun = mean(abun)) %>%
    data.frame %>% arrange(-mean.abun)
df.ab$phylum %>% unique()
phy = c(df.ab %>% filter(phylum != 'Remainder') %>% arrange(-mean.abun) %>% .$phylum %>% unique(), 'Remainder')

options(
  ggplot2.discrete.colour = pals::polychrome(25) %>% data.frame %>% .[3:25,1],
  ggplot2.discrete.fill = c("#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", "#B00068", "#AA0DFE", "#90AD1C", "#2ED9FF", "#DEA0FD",'#353636')
)

df.ab %>% mutate(phylum = factor(phylum, levels = phy)) %>%
 ggplot(aes(x = site, y = mean.abun, fill = phylum)) +
        geom_bar(stat = 'identity', color = '#353636', alpha = 0.8) + 
        facet_grid(.~sym) + 
        ylab("Relative abundance") +
        guides(fill=guide_legend(title="phylum")) + 
        scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.1))) + theme_bw()
ggsave(file = 'figures/relative-abundance_bacteria_site-sym_2017-Bulk_phylum.svg', width = 6, height = 3, dpi = 900, units = 'in')

########################################
########################################
df.ab = ps %>% 
    speedyseq::filter_tax_table(!is.na(phylum)) %>%
    speedyseq::tax_glom(taxrank = "phylum") %>%
    speedyseq::transform_sample_counts(function(x) x / sum(x)) %>%
    speedyseq::psmelt() %>% filter(date == '2022' & soil == 'Rhizosphere') %>%
    data.frame %>% 
    group_by(site, sym, phylum) %>% 
    mutate(median =  median(Abundance)) %>%
    mutate(phylum = ifelse(median <= 0.01, 'Remainder', phylum)) %>%
    data.frame %>%
    group_by(site, sym, phylum, Sample)%>%
    dplyr::summarize(abun=sum(Abundance))%>%
    group_by(site, sym, phylum)%>%
    dplyr::summarize(mean.abun = mean(abun)) %>%
    data.frame %>% arrange(-mean.abun)

df.ab$phylum %>% unique()
phy = c(df.ab %>% filter(phylum != 'Remainder') %>% arrange(-mean.abun) %>% .$phylum %>% unique(), 'Remainder')

options(
  ggplot2.discrete.colour = pals::polychrome(25) %>% data.frame %>% .[3:25,1],
  ggplot2.discrete.fill = c("#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", 
                            "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", '#353636')
)

df.ab %>% mutate(phylum = factor(phylum, levels = phy)) %>%
 ggplot(aes(x = site, y = mean.abun, fill = phylum)) +
        geom_bar(stat = 'identity', color = '#353636', alpha = 0.8) + 
        facet_grid(.~sym) + 
        ylab("Relative abundance") +
        guides(fill=guide_legend(title="phylum")) + 
        scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.1)))
ggsave(file = 'figures/relative-abundance_bacteria_site-sym_2022-Rhizosphere_phylum.svg', width = 6, height = 3, dpi = 900, units = 'in')

########################################
df.ab = ps %>% 
    speedyseq::filter_tax_table(!is.na(phylum)) %>%
    speedyseq::tax_glom(taxrank = "phylum") %>%
    speedyseq::transform_sample_counts(function(x) x / sum(x)) %>%
    speedyseq::psmelt() %>% filter(date == '2022' & soil == 'Bulk') %>%
    data.frame %>% 
    group_by(site, sym, phylum) %>% 
    mutate(median =  median(Abundance)) %>%
    mutate(phylum = ifelse(median <= 0.01, 'Remainder', phylum)) %>%
    data.frame %>%
    group_by(site, sym, phylum, Sample)%>%
    dplyr::summarize(abun=sum(Abundance))%>%
    group_by(site, sym, phylum)%>%
    dplyr::summarize(mean.abun = mean(abun)) %>%
    data.frame %>% arrange(-mean.abun)
df.ab$phylum %>% unique()
phy = c(df.ab %>% filter(phylum != 'Remainder') %>% arrange(-mean.abun) %>% .$phylum %>% unique(), 'Remainder')

options(
  ggplot2.discrete.colour = pals::polychrome(25) %>% data.frame %>% .[3:25,1],
  ggplot2.discrete.fill = c("#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", 
                            "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", '#353636')
)

df.ab %>% mutate(phylum = factor(phylum, levels = phy)) %>%
 ggplot(aes(x = site, y = mean.abun, fill = phylum)) +
        geom_bar(stat = 'identity', color = '#353636', alpha = 0.8) + 
        facet_grid(.~sym) + 
        ylab("Relative abundance") +
        guides(fill=guide_legend(title="phylum")) + 
        scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.1))) + theme_bw()
ggsave(file = 'figures/relative-abundance_bacteria_site-sym_2022-Bulk_phylum.svg', width = 6, height = 3, dpi = 900, units = 'in')
