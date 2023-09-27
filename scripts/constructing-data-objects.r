# Load packages
library(tidyverse)
library(speedyseq) # this is a package from github, the package "phyloseq" can be used instead. "Speedyseq" is a faster version of "phyloseq"
library(vegan)

# Load data
load("data/Brzostek_CABBI-seqtabnochim.RData")
load("data/Brzostek_CABBI-tax-table-object.RData")
load("data/Brzostek_CABBI-tree-object.RData")
df = read.csv(file = "data/mycorrhizal_kangi_brzostek_metadata.csv", header = TRUE, sep = ",") %>% 
        select(sample_id, date, site, sym, soil, myco) %>% 
        mutate(date = as.character(date))
# Calculate alpha diversity
df.ens = vegan::diversity(seqtab.nochim, index = "invsimpson") %>% data.frame('ens' = `.`) %>% mutate('sample_id' = rownames(.))
df.rich = vegan::estimateR(seqtab.nochim) %>% t() %>% data.frame %>% select(S.obs) %>% mutate('sample_id' = rownames(.), 'rich' = S.obs) %>% select(sample_id, rich)
df.alpha = left_join(df.ens, df.rich, by = "sample_id")

df.meta = left_join(df, df.alpha, by = "sample_id")
rownames(df.meta) = df.meta$sample_id

# Build phyloseq object
ps = phyloseq::phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
                        sample_data(df.meta), 
                        tax_table(tax.table), 
                        phy_tree(tree))
ps = prune_samples(sample_sums(ps) > 4000, ps)
df.meta = df.meta %>% filter(sample_id %in% sample_names(ps))
# Save data objects
save(df.meta, ps, file = "data/fernow-analysis-data.RData")
