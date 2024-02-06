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

permanova(ps %>% filter_sample_data(sym == "AM" & site == 'ws3'), "date * soil")
# sig diff between dates
plot.nmds(ps %>% filter_sample_data(sym == "AM" & site == 'ws3'), "date", 'date')


permanova(ps %>% filter_sample_data(sym == "AM" & site == 'ws7'), "date * soil")
# sig diff between dates
plot.nmds(ps %>% filter_sample_data(sym == "AM" & site == 'ws7'), "date", 'date')


permanova(ps %>% filter_sample_data(sym == "ECM" & site == 'ws3'), "date * soil")
# sig diff between dates
plot.nmds(ps %>% filter_sample_data(sym == "ECM" & site == 'ws3'), "date", 'date')


permanova(ps %>% filter_sample_data(sym == "ECM" & site == 'ws7'), "date * soil")
# NO sig diff between dates



ps %>% speedyseq::filter_sample_data(date == "2017") %>%
    speedyseq::filter_sample_data(sym == "AM" & site == "ws3" & soil == "Rhizosphere") %>%
    sample_names()




pwd.time = function(sym, site, soil){
    dist = speedyseq::otu_table(ps) %>% vegan::vegdist(., method = "bray") %>% as.matrix()
    
    n17 = ps %>% speedyseq::filter_sample_data(date == "2017") %>%
        speedyseq::filter_sample_data(sym == {{ sym }} & site == {{ site }} & soil == {{ soil }}) %>%
        sample_names()
    n22 = ps %>% speedyseq::filter_sample_data(date == "2022") %>%
        speedyseq::filter_sample_data(sym == {{ sym }} & site == {{ site }} & soil == {{ soil }}) %>%
        sample_names()

    dfpwd = dist[n17, n22] %>% 
        rowMeans() %>% 
        data.frame() %>% 
        magrittr::set_colnames("pwd") %>%
        #summarize(mean = mean(`.`), var = var(`.`), se = sd(`.`)/sqrt(n())) %>%
        mutate(sym = sym, site = site, soil = soil) %>% 
        select(pwd, sym, site, soil)
    
    return(dfpwd)
}

df.pw = data.frame(expand_grid(sym = c("AM", "ECM"), site = c("ws3", "ws7"), soil = c("Rhizosphere", "Bulk"))) %>%
    plyr::mdply(., pwd.time) %>% bind_rows()

lm(pwd ~ sym * site * soil, data = df.pw) %>% 
    emmeans::emmeans(., ~ sym * site * soil) %>%
    multcomp::cld(., Letters = letters) %>%
    data.frame() %>%
    ggplot(aes(x = site, y = emmean, fill = sym)) +
        geom_linerange(aes(x = site, ymin = emmean - SE, ymax = emmean + SE), 
                            position = position_dodge(width = 0.6)) +
        geom_point(shape = 21, size = 4, position = position_dodge(width = 0.6)) +
        geom_text(aes(label = .group, y = ((emmean + SE) + ((emmean + SE) * 0.05))), position = position_dodge(width = 0.6), vjust = -1) +
        facet_grid(.~soil)


# Samples are not different between 2017 and 2022 at a compositional level.

bind_rows(pwd.time("AM", "ws3", "Rhizosphere"),
pwd.time("AM", "ws7", "Rhizosphere"),
pwd.time("ECM", "ws3", "Rhizosphere"),
pwd.time("ECM", "ws7", "Rhizosphere"),
pwd.time("AM", "ws3", "Bulk"),
pwd.time("AM", "ws7", "Bulk"),
pwd.time("ECM", "ws3", "Bulk"),
pwd.time("ECM", "ws7", "Bulk")) %>%
    mutate(sym = factor(sym, levels = c("AM", "ECM")),
           site = factor(site, levels = c("ws3", "ws7")),
           soil = factor(soil, levels = c("Rhizosphere", "Bulk"))) %>%
    ggplot(aes(x = site, y = mean, fill = sym)) +
        geom_linerange(aes(x = site, ymin = mean - se, ymax = mean + se), 
                            position = position_dodge(width = 0.6)) +
        geom_point(shape = 21, size = 4, position = position_dodge(width = 0.6)) +
        facet_grid(.~soil)


