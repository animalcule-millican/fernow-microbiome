library(tidyverse)
library(vegan)
library(speedyseq)
load("data/fernow-analysis-data.RData")
##################################################
######## Pairwise distances and variances ########
##################################################
pairwise.var = function(PS, site, sym, soil, date){
    dist = PS %>% 
            speedyseq::otu_table() %>% 
            vegdist(., method = "bray") %>% 
            as.matrix
    
    sam.names = PS %>% filter_sample_data(site == {{ site }} & 
                                            sym == {{ sym }} & 
                                            soil == {{ soil }} & 
                                            date == {{ date }}) %>%
                    phyloseq::sample_names()

    df.var = dist[sam.names, sam.names] %>% 
        matrixStats::rowVars() %>% 
        data.frame('var' = `.`) %>% 
        mutate(sample_id = rownames(.)) %>%
        mutate(site = {{ site }},
                sym = {{ sym }},
            soil = {{ soil }},
            date = {{ date }})
    return(df.var)
}

df.var = data.frame(expand_grid('site' = df.meta$site %>% unique(), 
                    'sym' = df.meta$sym %>% unique(), 
                    'soil' = df.meta$soil %>% unique(), 
                    'date' = df.meta$date %>% unique())) %>% 
    plyr::mdply(., pairwise.var, PS = ps)

df.var %>%
    lm(var ~ site * sym * soil * date, data = .) %>%
    car::Anova()

df.stat = df.var %>%
    lm(var ~ site * sym * soil * date, data = .) %>%
    emmeans::emmeans(., ~ date + sym | soil + site) %>%
    multcomp::cld(., Letters = letters) %>%
    data.frame() %>%
    mutate(group = toupper(trimws(`.group`))) %>%
    select(-.group)
##################################################
df.var %>%
    ggplot(aes(x = sym, y = 1-var, color = date)) +
    geom_point(position = position_dodge(width = 0.6)) +
    geom_text(data = df.stat, aes(x = sym, y = ((1-upper.CL)+((1-upper.CL) * 0.05)), label = group),
                position = position_dodge(width = 0.6), show.legend = FALSE) +
    facet_grid(site ~ soil) +
    theme_bw() + ylab("Community Stability Index (1 - pwd variance)") +
    theme(axis.title.x = element_blank())
ggsave("figures/community-stability_Sym-Site.svg", width = 3, height = 3, units = "in", dpi = 600)

df.var %>%
    ggplot(aes(x = soil, y = 1-var, color = sym)) +
    geom_point(position = position_dodge(width = 0.6)) +
    geom_text(data = df.stat, aes(x = soil, y = ((1-upper.CL)+((1-upper.CL) * 0.05)), label = group),
                position = position_dodge(width = 0.6), show.legend = FALSE) +
    facet_grid(site ~ date) +
    theme_bw() + ylab("Community Stability Index (1 - pwd variance)") +
    theme(axis.title.x = element_blank())
ggsave("figures/community-stability_Sym-Soil.svg", width = 3, height = 3, units = "in", dpi = 600)

df.var %>%
    ggplot(aes(x = site, y = 1-var, color = sym)) +
    geom_point(position = position_dodge(width = 0.6)) +
    geom_text(data = df.stat, aes(x = site, y = ((1-upper.CL)+((1-upper.CL) * 0.05)), label = group),
                position = position_dodge(width = 0.6), show.legend = FALSE) +
    facet_grid(soil ~ date) +
    theme_bw() + ylab("Community Stability Index (1 - pwd variance)") +
    theme(axis.title.x = element_blank())
ggsave("figures/community-stability_Sym-Site.svg", width = 3, height = 3, units = "in", dpi = 600)
