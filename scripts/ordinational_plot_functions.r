project_theme <- function(base_size = 12, font = 'Times') {
  theme_bw(base_size = base_size, base_family = font) %+replace%
    theme(
      # Plot title frmating
      text = element_text(family = font, face = 'plain', size = rel(1)),
      plot.title = element_text(margin = margin(0,0,5,0), hjust = 0),
      # Grid and border
      axis.title = element_text(size = rel(1.0)),
      panel.grid.major = element_line(color = 'lightgray', linewidth = 0.2),
      panel.grid.minor = element_blank(),
      # Axis formatting
      axis.line = element_line(color = "black"),
      # Legend formatting
      #legend.position="none",
      # Facet formatting
      strip.background = element_rect(fill="transparent", color = "transparent"),
      strip.text = element_text(colour = 'black')
    )
}

options(
  ggplot2.discrete.colour = c("#F6222E", "#3283FE"),
  ggplot2.discrete.fill = c("#F6222E", "#3283FE")
)

theme_set(project_theme())

permanova = function(PS, var){
    fml = formula(paste0("asv ~ ", var))
    asv = PS %>% otu_table %>% vegan::vegdist(., method = "bray")
    sam = PS %>% sample_data() %>% data.frame()
    vegan::adonis2(fml, data = sam, permutations = 999)
}

gordiplot <- function(ord, sample.df, group, shape = c('ellipse', 'hull', 'spider'), kind = c("sd", "se", "ehull")) {
  require('tidyverse')
  require('vegan')
  conf=NULL
  show.group="all"
  scaling = 1
  # Get site coordinates to plot.
  df_ord <- vegan::scores(ord, display = "sites", scaling=scaling) %>%
  data.frame() %>% mutate(sample_id = rownames(.)) %>% 
  left_join(., sample.df %>% select(sample_id, {{ group }})) %>%
  rename('Group' = {{ group }}) %>% 
  select(NMDS1, NMDS2, Group, sample_id) %>%
  magrittr::set_colnames(c('x','y','Group','sample_id'))

  group = sample.df %>% select(all_of(group)) %>% .[[group]]
  group <- as.factor(group)
  if (show.group[1]=="all") {
    show.group <- as.vector(levels(group))
  }
  # Add ellipses.
  if (shape == 'ellipse') {
    # Get ellipse centers to annotate.
  df_mean.ord <- aggregate(df_ord[,1:2], by=list(df_ord$Group),mean)
  colnames(df_mean.ord) <- c("Group", "x", "y")
  df_mean.ord <- df_mean.ord[df_mean.ord$Group %in% show.group, ]

  # Get parameters from the ordiellipse function.
  if (is.null(conf)) {
    rslt <- vegan::ordiellipse(ord, group=group, display = "sites", scaling=scaling, kind = kind, show.group = show.group, draw = "none")
  } else {
      rslt <- vegan::ordiellipse(ord, group=group, display = "sites", scaling=scaling, kind = kind, show.group = show.group, draw = "none", conf = conf)
  }

  # Get points to plot for the ellipses.
  df_ellipse <- data.frame()
  for(g in show.group) {
    df_ellipse <- rbind(df_ellipse, cbind(as.data.frame(with(df_ord[df_ord$Group==g,],
    vegan:::veganCovEllipse(rslt[[g]]$cov,rslt[[g]]$center, rslt[[g]]$scale))),Group=g))
  }
  colnames(df_ellipse) <- c("x", "y", "Group")
  df_ellipse <- df_ellipse[ , c(3,1,2)]
    return(df_ellipse)
  }


  # Add hulls.
  if (shape == 'hull') {
    # Make data frame for hulls.
  rslt.hull <- vegan::ordihull(ord, group = group, scaling = scaling, show.group = show.group, draw = "none")
  df_hull <- data.frame()
  df_temp <- data.frame()
  for (g in show.group) {
    x <- rslt.hull[[g]][ , 1]
    y <- rslt.hull[[g]][ , 2]
    Group <- rep(g, length(x))
    df_temp <- data.frame(Group = Group, x=x, y=y)
    df_hull <- rbind(df_hull, df_temp)
  }
    return(df_hull)
  
  }

  # Add spiders.
  
  if (shape == 'spiders') {
    # Make a data frame for the spiders.
  df_spiders <- df_ord
  df_spiders$cntr.x <- NA
  df_spiders$cntr.y <- NA
  for (g in show.group) {
    df_spiders[which(df_spiders$Group==g), 4:5] <- df_mean.ord[which(df_mean.ord==g), 2:3]
  }
  df_spiders <- df_spiders[ , c(3,4,5,1,2)]
  df_spiders <- df_spiders[order(df_spiders$Group), ]
  df_spiders <- df_spiders[df_spiders$Group %in% show.group, ]
    return(df_spiders)
    
  }

}
beta.base.plot = function(df,df2,df3, color){
    require(tidyverse)
    require(vegan)
    require(phyloseq)
    dat = df['points'] %>% data.frame %>% mutate(sample_id = row.names(.)) %>%
    left_join(., df2 %>% select(sample_id, {{ color }}), by = 'sample_id') 
    plot = dat %>% 
    ggplot(aes(x = points.MDS1, y = points.MDS2, color = {{ color }})) +
	  geom_path(data = df3,
		    aes(x = x, y = y, group = Group, color = Group), 
		    linetype = 6, alpha = 0.8, linewidth = 0.5, show.legend = FALSE) +
	  geom_point(aes(fill = {{ color }}), size = 2, shape = 21, color = 'black') +
	  xlab('nMDS1') +
	  ylab('nMDS2')
    return(plot)
}
plot.nmds = function(PS, var, color_var){
    require(tidyverse)
    require(phyloseq)
    fml = formula(paste0('comm ~ ', var))
    ps.tf = PS %>%
        phyloseq::transform_sample_counts(., function(x) sqrt(x/sum(x)))
    sam = phyloseq::sample_data(ps.tf) %>% data.frame
    asv = ps.tf %>% phyloseq::transform_sample_counts(., function(x) sqrt(x/sum(x))) %>% phyloseq::otu_table(.)
    comm = vegan::vegdist(asv, method = 'bray')
    nmds = vegan::metaMDS(comm, distance = 'bray', maxit = 9999, k = 2, trymax = 10000, try = 50)
    hull = gordiplot(nmds, sam, {{ color_var }}, shape = 'hull', kind = 'ehull')
    fg = beta.base.plot(nmds, sam, hull, !!sym(color_var))
    print(adonis2(fml, data = sam, permutations = 9999))
    return(fg)
}
