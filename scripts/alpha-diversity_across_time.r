library("tidyverse")
head(df.meta)
unique(df.meta$date)
unique(df.meta$site)
unique(df.meta$sym)

# AMF dominant samples

## Watershed 3
df.meta %>% filter(sym == "AM" & site == 'ws3') %>% 
    lm(rich ~ soil * date, data = .) %>%
    car::Anova()
df.meta %>% filter(sym == "AM" & site == 'ws3') %>% 
    lm(ens ~ soil * date, data = .) %>%
    car::Anova()
### No significant differences in AMF dominated sample richness or diversity between 2017 and 2022

## Watershed 7
df.meta %>% filter(sym == "AM" & site == 'ws7') %>% 
    lm(rich ~ soil * date, data = .) %>%
    car::Anova()
df.meta %>% filter(sym == "AM" & site == 'ws7') %>% 
    lm(ens ~ soil * date, data = .) %>%
    car::Anova()
### Significant difference in richness between 2017 and 2022 for AMF samples in ws7
#### This is not the case for ens diversity



# ECM dominant samples

## Watershed 3
df.meta %>% filter(sym == "ECM" & site == 'ws3') %>% 
    lm(rich ~ soil * date, data = .) %>%
    car::Anova()
df.meta %>% filter(sym == "ECM" & site == 'ws3') %>% 
    lm(ens ~ soil * date, data = .) %>%
    car::Anova()
### There is a significant difference in richness between 2017 and 2022 for ECM samples in ws3
#### This is not the case for ens diversity

## Watershed 7
df.meta %>% filter(sym == "ECM" & site == 'ws7') %>% 
    lm(rich ~ soil * date, data = .) %>%
    car::Anova()
df.meta %>% filter(sym == "ECM" & site == 'ws7') %>% 
    lm(ens ~ soil * date, data = .) %>%
    car::Anova()
### No observed differences between 2017 and 2022 for ECM samples in ws7