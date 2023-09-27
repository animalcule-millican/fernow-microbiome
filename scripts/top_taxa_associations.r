# load packages
library(tidyverse)
library(speedyseq)
# load data
load("data/fernow-analysis-data.RData")
#######################################
# Identifying top taxa related to differences in communitiy composition and mycorrhizal root colonization.


get_coef = function(PS){
    sam = sample_data(PS) %>% data.frame() %>% na.omit()
    
    rel_abund <- PS %>% transform_sample_counts(., function(x) x / sum(x)) %>% filter_sample_data(sample_id %in% sam$sample_id) %>%
                    otu_table
    
    tax.df = tax_table(PS) %>% data.frame %>% mutate(seqid = rownames(.))
    
    permanova <- vegan::adonis(rel_abund ~ myco, data = sample_data(PS) %>% data.frame())
    
    coef <- coefficients(permanova) %>% data.frame
    
    top.tax = coef %>% 
    t() %>% 
    data.frame %>% 
    filter(myco != 0) %>%
    arrange(-(abs(myco))) %>%
    mutate(asv = rownames(.), coef = myco) %>% 
    head(n = 20) %>%
    select(asv, coef) %>% 
    data.frame() %>%
    arrange(-coef) %>% 
    left_join(., tax.df, by = c('asv' = 'seqid')) %>%
    #mutate(asv = seqid) %>%
    magrittr::set_rownames(NULL)
    return(top.tax)
}


df.co1 = get_coef(ps %>% filter_sample_data(date == '2017' & site == "7_ecosystem_recovery"))
df.co2 = get_coef(ps %>% filter_sample_data(date == '2017' & site != "7_ecosystem_recovery"))
df.co3 = get_coef(ps %>% filter_sample_data(date == '2022' & site == "7_ecosystem_recovery"))
df.co4 = get_coef(ps %>% filter_sample_data(date == '2022' & site != "7_ecosystem_recovery"))
df.co4$genus %>% unique()

df.gen = ps %>% tax_glom(taxrank = "genus") %>% transform_sample_counts(., function(x) x / sum(x)) %>% speedyseq::psmelt()

df.gen %>% 
    filter(genus == 'Actinoallomurus' & date == '2017') %>%
    ggplot(aes(x = myco, y = Abundance, color = sym)) +
    geom_point() +
    stat_smooth(method = 'lm', se = FALSE) +
    ggpubr::stat_cor()

df.gen %>% 
    filter(genus == 'Aquisphaera' & date == '2017') %>%
    ggplot(aes(x = myco, y = Abundance, color = sym)) +
    geom_point() +
    stat_smooth(method = 'lm', se = FALSE) +
    ggpubr::stat_cor()
df.gen %>% 
    filter(genus == 'Mycobacterium' & date == '2017') %>%
    ggplot(aes(x = myco, y = Abundance, color = sym)) +
    geom_point() +
    stat_smooth(method = 'lm', se = FALSE) +
    ggpubr::stat_cor()

df.gen %>% 
    filter(genus == 'Diplorickettsia' & date == '2017') %>%
    ggplot(aes(x = myco, y = Abundance, color = sym)) +
    geom_point() +
    stat_smooth(method = 'lm', se = FALSE) +
    ggpubr::stat_cor()

df.gen %>% 
    filter(genus == 'Rhodoplanes' & date == '2017') %>%
    ggplot(aes(x = myco, y = Abundance, color = sym)) +
    geom_point() +
    stat_smooth(method = 'lm', se = FALSE) +
    ggpubr::stat_cor()

df.gen %>% 
    filter(genus == 'Gemmata' & date == '2017') %>%
    ggplot(aes(x = myco, y = Abundance, color = sym)) +
    geom_point() +
    stat_smooth(method = 'lm', se = FALSE) +
    ggpubr::stat_cor()

df.gen %>% 
    filter(genus == 'Methylocapsa' & date == '2017') %>%
    ggplot(aes(x = myco, y = Abundance, color = sym)) +
    geom_point() +
    stat_smooth(method = 'lm', se = FALSE) +
    ggpubr::stat_cor()
# Interesting, Methylocapsa acidophila (an N-fixing bacteria) is found in wood colonized by Basidiomycete Hypholoma fasciculare. 
# This assoication is thought to lead to CH4 production by the fungi (citing references: DOI:10.1111/j.1574-6941.2007.00425.x)
#  https://doi.org/10.1038/ncomms2049
## Look for links between AMF and methanotrophs in literature.