---
Author: Michael D Millican
---

```{r, echo=FALSE include=FALSE}
library(tidyverse)
library(speedyseq)
# load data
load("data/fernow-analysis-data.RData")
```
# How is the alpha diversity of bacterial communities impacted by fertilzation history?

"""
Is there a difference in bacterial richness between soil types, within each watershed, year, and symbiotroph dominance?
"""
df.meta %>% 
    lm(rich ~ site * sym * soil * date, data = .) %>%
    emmeans::emmeans(., ~ soil | date * sym * site) %>%
    multcomp::cld(., Letters = letters)

"""
No. 
There is no significant difference in bacterial richness between rhizophere and bulk soils, within each watershed, year, or symbiotroph ecosystem.
"""

"""
Is there a difference in bacterial diversity between soil types, within each watershed, year, and symbiotroph dominance?
"""
df.meta %>% 
    lm(ens ~ site * sym * soil * date, data = .) %>%
    emmeans::emmeans(., ~ soil | date * sym * site) %>%
    multcomp::cld(., Letters = letters)

"""
No. 
There is no significant difference in bacterial diversity between rhizophere and bulk soils, within each watershed, year, or symbiotroph ecosystem.
"""

# Conclusion:
    - The soil type does not play a significant role in the determination of bacterial alpha diversity at any watershed site, sampling year, or symbiotroph dominated ecosystem.


"""
Is there a difference in alpha diversity between sampling years, within each watershed, soil type, and symbiotroph dominance?
"""
df.meta %>% 
    lm(rich ~ site * sym * soil * date, data = .) %>%
    emmeans::emmeans(., ~ date | soil * sym * site) %>%
    multcomp::cld(., Letters = letters)
- ECM
    - ws3
        - Bulk
            - no dif
        - Rhizosphere
            - sig dif
                - 2017 > 2022
    - ws7
        - Bulk
            - no difference
        - Rhizosphere
            - no difference
- AM
    - ws3
        - Bulk
            - no dif
        - Rhizosphere
            - no dif
    - ws7
        - Bulk
            - sig dif
                - 2017 < 2022
        - Rhizosphere
            - sig dif 
                - 2017 > 2022


df.meta %>% 
    lm(ens ~ site * sym * soil * date, data = .) %>%
    emmeans::emmeans(., ~ date | soil * sym * site) %>%
    multcomp::cld(., Letters = letters)

- No differences observed in bacterial ENSpie diversity.


###

df.meta %>% 
    lm(rich ~ site * sym * soil * date, data = .) %>%
    emmeans::emmeans(., ~ site * date * soil * sym) %>%
    multcomp::cld(., Letters = letters)
