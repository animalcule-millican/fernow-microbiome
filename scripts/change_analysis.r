library(tidyverse)
library(speedyseq)
# load data
load("data/fernow-analysis-data.RData")

# Changes in alpha diversity between 2017 and 2022
df = read.csv("data/mycorrhizal_kangi_brzostek_metadata_extra_data.csv", header = TRUE) %>%
        select(sample_id, sample_label)

split_ids = function(x){
    strsplit(x, "(?<=AM)|(?<=ECM)", perl = TRUE)[[1]][2] %>%
    strsplit(., " ", perl = TRUE) %>%
    unlist() %>%
    .[1]

}

df$rep = df$sample_label %>% lapply(., split_ids) %>% unlist()

df.meta = df.meta %>%
    left_join(., df, by = "sample_id")
df.17 = df.meta %>% filter(date == "2017") %>%
            arrange(date, site, soil, sym, rep)
df.22 = df.meta %>% filter(date == "2022") %>%
            arrange(date, site, soil, sym, rep)
df.22 %>% select(sample_id, site, sym, soil, rep) %>%
    mutate(rich.chng = df.22$rich - df.17$rich,
           ens.chng = df.22$ens - df.17$ens)

df.meta %>% group_by(site, sym, soil, rep) %>%
    summarise(rich.chg = rich[date == '2022'] - rich[date == '2017'],
                ens.chg = ens[date == '2022'] - ens[date == '2017']) %>%
    lm(rich.chg ~ site * sym * soil, data = .) %>%
    car::Anova()



df.meta %>% group_by(site, sym, soil, rep) %>%
    summarise(rich.chg = rich[date == '2022'] - rich[date == '2017'],
                ens.chg = ens[date == '2022'] - ens[date == '2017']) %>%
    lm(ens.chg ~ site * sym * soil, data = .) %>%
    car::Anova()
