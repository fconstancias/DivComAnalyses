#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'library(phyloseq)
#'data(esophagus)
#'phyloseq_compute_bdiv(esophagus, 100) -> dist
#'
#'
#'
#'
#'
#'
#'###########################################################################################

phyloseq_compute_bdiv <- function(phylo_tmp, 
                                  norm, 
                                  phylo,
                                  seed)
{
  require(tidyverse)
  require(ape)
  require(phyloseq)
  require(abdiv)
  require(GUniFrac)
  
  set.seed(seed)
  if(norm == "pc")
  {
    phylo_tmp %>% 
      transform_sample_counts(function(x) {x/sum(x)} * 100) -> phylo_tmp
  } else{
    phylo_tmp = phylo_tmp
  }
  
  if (phylo == TRUE)
  {
    
    #https://github.com/joey711/phyloseq/issues/936
    phy_tree(phylo_tmp) <- ape::multi2di(phy_tree(phylo_tmp))
    
    dist_methods <- c("bray","bjaccard", "wjaccard", "uunifrac", "wunifrac"); dist_unif <- c("d_0", "d_0.5")
    
    dlist <- vector("list", length(dist_methods) + length(dist_unif))
    names(dlist) = c(dist_methods, dist_unif)
    
    # compute phyloseq::distance distances
    for( i in dist_methods ){
      if ( i == "bjaccard"){ # binary disances computed using vegan requires to specify binary = TRUE
        phylo_tmp %>% 
          phyloseq::distance(method = "jaccard", 
                             binary = TRUE) ->  dlist[[i]]
      } else if  ( i == "wjaccard"){ # binary disances computed using vegan requires to specify binary = TRUE
        phylo_tmp %>% 
          phyloseq::distance(method = "jaccard", 
                             binary = FALSE) ->  dlist[[i]]
      } else{
        phylo_tmp %>% 
          phyloseq::distance(method = i) ->  dlist[[i]]
      }
    }
    set.seed(seed)
    unifracs <- GUniFrac::GUniFrac(phylo_tmp %>% otu_table() 
                                   %>% t(),
                                   phy_tree(phylo_tmp), 
                                   alpha=c(0, 0.5))$unifracs
    
    for( i in dist_unif ){
      
      dlist[[i]] <- unifracs[, , i] %>% as.dist()
    }
  }else{
    dist_methods <- c("bray", "sorensen","bjaccard", "wjaccard")
    
    dlist <- vector("list", length(dist_methods))
    names(dlist) = c(dist_methods)
    
    # compute phyloseq::distance distances
    for( i in dist_methods ){
      if ( i== "bjaccard"){ # binary disances computed using vegan requires to specify binary = TRUE
        phylo_tmp %>% 
          phyloseq::distance(method = "jaccard", 
                             binary = TRUE) ->  dlist[[i]]
      } else if  ( i == "wjaccard"){ # binary disances computed using vegan requires to specify binary = TRUE
        phylo_tmp %>% 
          phyloseq::distance(method = "jaccard", 
                             binary = FALSE) ->  dlist[[i]]
      } else if  ( i == "sorensen"){ # binary disances computed using vegan requires to specify binary = TRUE
        phylo_tmp %>% 
          phyloseq::distance(method = "bray", 
                             binary = TRUE) ->  dlist[[i]]
      }else{
        phylo_tmp %>% 
          phyloseq::distance(method = i) ->  dlist[[i]]
      }
    }
    
  }
  return(dlist)
  
}

###########################################################################################
phyloseq_plot_bdiv <- function(ps_rare,
                               dlist,
                               m = "PCoA",
                               seed,
                               axis1 = axis1,
                               axis2 = axis2)
{
  if (m == "CoDa")
  {
    
    d.czm <- zCompositions::cmultRepl(t(otu_table(ps_rare)),  label=0, method="CZM")
    d.clr <- t(apply(d.czm, 1, function(x){log(x) - mean(log(x))})) %>% as.matrix() 
    
    ps_rare@otu_table = NULL
    
    ps_rare@otu_table = otu_table((d.clr), taxa_are_rows = FALSE)
    
    # ps_clr <- microbiome::transform(ps_rare, "clr")
    ord_clr <- phyloseq::ordinate(ps_rare, "RDA")
    
    #Scale axes and plot ordination
    phyloseq::plot_ordination(ps_rare, ord_clr, type="samples") -> p
      # coord_fixed(ord_clr$CA$eig[2] / sum(ord_clr$CA$eig) / ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)) 
    #https://github.com/ggloor/CoDa_microbiome_tutorial/wiki/Part-1:-Exploratory-Compositional-PCA-biplot
    
    
    out <- list("PCA" = p, 
                "physeq_clr" = ps_rare)
    
    return(out)
      }else{
    
  plot_list <- vector("list", length(dlist))
  names(plot_list) =  names(dlist)
  
  for( i in dlist %>% names){
    print(i)
    set.seed(seed)
    
    if(m == "TSNE")
    {
      as.matrix(dlist[[i]])[sample_names(ps_rare),sample_names(ps_rare)] %>%
        as.dist() -> dlist[[i]]
      
      # https://microbiome.github.io/microbiome/Ordination.html
      # Run TSNE
      tsne_out <- Rtsne::Rtsne(dlist[[i]], dims = 2, perplexity = 5, verbose = T) 
      proj <- tsne_out$Y %>% data.frame()
      
      rownames(proj) <- rownames(t(otu_table(ps_rare)))
      
      proj2 <- cbind(proj, sample_data(ps_rare))
      
      # rownames(proj) == sample_data(ps_rare)$SampleID
      
      # microbiome::plot_landscape(proj, legend = T, size = 1) 
      p <- phyloseq::plot_ordination(ps_rare, 
                                     proj) + geom_point(data = proj2, aes_string(colour = "Oral_Site", shape = "Health_status"), size = 4) + theme_bw() #+
      # scale_shape_manual(values=seq(0,15))
      plot_list[[i]] = p
      
    }else{
      
      as.matrix(dlist[[i]])[sample_names(ps_rare),sample_names(ps_rare)] %>%
        as.dist() -> dlist[[i]]
      
      # Calculate ordination
      iMDS  <- ordinate(ps_rare, 
                        m,
                        distance = dlist[[i]])
      
      # Create plot, store as temp variable, p
      p <- phyloseq::plot_ordination(ps_rare, iMDS, 
                                     axes = axis1:axis2)
      # Add title to each plot
      if(m == "NMDS")
      {
        p <- p + ggtitle(paste0(m," using distance method ", i, "\n",
                                " NMDS 2d stress = ", iMDS$grstress %>% round(2))) +
          geom_point(size = 4) + theme_bw() #+  
        # ggrepel::geom_text_repel(cex=2.5,aes(label=sample))
        plot_list[[i]] = p
      }
      if(m == "PCoA")
      {
        p <- p + ggtitle(paste0(m," using distance method ", i)) +
          geom_point(size = 4) + theme_bw() #+  
        # ggrepel::geom_text_repel(cex=2.5,aes(label=sample))}
        # Save the graphic to file.
        plot_list[[i]] = p
        # # Save the pairwise permanova
        # tmp_list_2[[i]] <- vegan::adonis(dlist[[i]] ~ get_variable(tmp, color))$aov.tab %>% data.frame()
        # # Save betadisper
        # tmp_list_3[[i]]  <- vegan::permutest(vegan::betadisper(dlist[[i]], get_variable(tmp, color)))$tab$`Pr(>F)`[1]
        # TW https://github.com/alekseyenko/WdStar/blob/master/16S_alone_taxa_of_interest.html
        # tmp_list_4[[i]]  <- Tw2.test(dlist[[i]], get_variable(tmp, color)) %>% as.data.frame() #%>% mutate(Day = day)
        
      }
    }
  }
  
  return(plot_list)
}
}
###########################################################################################

phyloseq_plot_PCoA_3d <- function(ps_rare,
                                  dlist,
                                  m = "PCoA",
                                  seed,
                                  color,
                                  shape)
{
  plot_list <- vector("list", length(dlist))
  names(plot_list) =  names(dlist)
  
  for( i in dlist %>% names){
    print(i)
    set.seed(seed)
    
    if ( m == "PCoA")
    {
      # Calculate ordination
      ord <- cmdscale(dlist[[i]], k = 3, eig =T) 
      ordata <- as.data.frame(ord$points)
      
      rownames(ordata) == sample_data(ps_rare)$SampleID
      ordata$Sample <- rownames(ordata)
      ordata <- cbind(ordata, sample_data(ps_rare))
      
      pty_pcoa <- plotly::plot_ly(ordata, x= ~V1, y=~V2, z = ~V3, 
                                  color = as.formula(paste0("~",paste0(color))), # #Description
                                  symbol = as.formula(paste0("~",paste0(shape))),
                                  colors = ggpubr::get_palette(palette = "npg", 
                                                               length(levels(get_variable(ps_rare,paste(color))))),
                                  symbols = c('circle','o'),
                                  marker = list(size = 10))  %>% 
        plotly::add_markers() %>%
        plotly::layout(title = paste0(m, ' - ', i),
                       scene = list(xaxis = list(title = paste0("PCoA 1 : ",round(ord$eig[1]/sum(ord$eig)*100,1),"%")),
                                    yaxis = list(title =  paste0("PCoA 2 : ",round(ord$eig[2]/sum(ord$eig)*100,1),"%")),
                                    zaxis = list(title =  paste0("PCoA 3 : ",round(ord$eig[3]/sum(ord$eig)*100,1),"%"))))
    }
    if ( m == "NMDS")
    {
      # Calculate ordination
      ord <- vegan::metaMDS(dlist[[i]], k = 3) 
      ordata <- as.data.frame(ord$points)
      
      rownames(ordata) == sample_data(ps_rare)$SampleID
      ordata$Sample <- rownames(ordata)
      ordata <- cbind(ordata, sample_data(ps_rare))
      
      pty_pcoa <- plotly::plot_ly(ordata, x= ~MDS1, y=~MDS2, z = ~MDS3, 
                                  color = as.formula(paste0("~",paste0(Group))), # #Description
                                  colors = ggpubr::get_palette(palette = "npg", 
                                                               length(levels(get_variable(ps_rare,paste(Group))))))  %>% 
        plotly::add_markers() %>%
        plotly::layout(title = paste0(m, ' - ', i),
                       scene = list(xaxis = list(title = paste0("NMDS 1")),
                                    yaxis = list(title =  paste0("NMDS 2")),
                                    zaxis = list(title =   paste0("NMDS 3"))))
    }
    plot_list[[i]] = pty_pcoa
  }
  return(plot_list)
}
###########################################################################################

calc_pairwise_permanovas_strata <- function(dm, metadata_map, compare_header, n_perm, strat) {
  require(mctoolsr)
  
  as.matrix(dm)[sample_names(physeq),sample_names(physeq)] %>%
    as.dist() -> dm
  
  comp_var = metadata_map[, compare_header]
  comp_pairs = combn(levels(comp_var), 2)
  pval = c()
  R2 = c()
  for (i in 1:ncol(comp_pairs)) {
    pair = comp_pairs[, i]
    dm_w_map = list(dm_loaded = dm, map_loaded = metadata_map)
    dm_w_map$map_loaded$in_pair = comp_var %in% pair
    dm_w_map_filt = mctoolsr::filter_dm(dm_w_map, filter_cat = "in_pair", 
                                        keep_vals = TRUE)
    
    if (strat %in% colnames(metadata_map)){
      
      if (!missing(n_perm)) {
        m = vegan::adonis(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[, 
                                                                             compare_header], permutations = n_perm,
                          strata = dm_w_map_filt$map_loaded[, 
                                                            strat])
      }
      else {
        m = vegan::adonis(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[, 
                                                                             compare_header],
                          strata = dm_w_map_filt$map_loaded[, 
                                                            strat])
      }
    }else{
      if (!missing(n_perm)) {
        m = vegan::adonis(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[, 
                                                                             compare_header], permutations = n_perm)
      }
      else {
        m = vegan::adonis(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[, 
                                                                             compare_header])
      }
    }
    pval = c(pval, m$aov.tab$`Pr(>F)`[1])
    R2 = c(R2, m$aov.tab$R2[1])
  }
  results = data.frame(t(comp_pairs), R2, pval)
  results$pvalBon = pval * length(pval)
  results$pvalFDR = round(pval * (length(pval)/rank(pval, ties.method = "average")), 
                          3)
  
  detach("package:mctoolsr", unload=TRUE)
  return(results)
}

###########################################################################################

physeq_pairwise_permanovas <- function(dm, physeq, compare_header, n_perm, strat) {
  require(mctoolsr)
  
  as.matrix(dm)[sample_names(physeq),sample_names(physeq)] %>%
    as.dist() -> dm
  
  physeq %>%
    sample_data() %>%
    data.frame() -> metadata_map
  
  comp_var = as.factor(metadata_map[, compare_header])
  comp_pairs = combn(levels(comp_var), 2)
  pval = c()
  R2 = c()
  for (i in 1:ncol(comp_pairs)) {
    pair = comp_pairs[, i]
    dm_w_map = list(dm_loaded = dm, map_loaded = metadata_map)
    dm_w_map$map_loaded$in_pair = comp_var %in% pair
    dm_w_map_filt = mctoolsr::filter_dm(dm_w_map, filter_cat = "in_pair", 
                                        keep_vals = TRUE)
    
    if (strat %in% colnames(metadata_map)){
      
      if (!missing(n_perm)) {
        m = vegan::adonis(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[, 
                                                                             compare_header], permutations = n_perm,
                          strata = dm_w_map_filt$map_loaded[, 
                                                            strat])
      }
      else {
        m = vegan::adonis(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[, 
                                                                             compare_header],
                          strata = dm_w_map_filt$map_loaded[, 
                                                            strat])
      }
    }else{
      if (!missing(n_perm)) {
        m = vegan::adonis(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[, 
                                                                             compare_header], permutations = n_perm)
      }
      else {
        m = vegan::adonis(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[, 
                                                                             compare_header])
      }
    }
    pval = c(pval, m$aov.tab$`Pr(>F)`[1])
    R2 = c(R2, m$aov.tab$R2[1])
  }
  results = data.frame(t(comp_pairs), R2, pval)
  results$pvalBon = pval * length(pval)
  results$pvalFDR = round(pval * (length(pval)/rank(pval, ties.method = "average")), 
                          3)
  
  detach("package:mctoolsr", unload=TRUE)
  
  return(results)
}

###########################################################################################

in_vitro_mIMT_STABvsTreat <- function(physeq,
                                      dlist,
                                      group,
                                      formula,
                                      strata,
                                      group_color,
                                      group_shape,
                                      group_alpha)
{
  
  ps <- prune_samples(get_variable(physeq, "treatment") == group,
                physeq)
  
  # ps <- physeq %>%
  #   subset_samples(treatment %in% group) 
  
  lapply(
    dlist,
    FUN = phyloseq_adonis,
    physeq = ps,
    # subset_samples(!treatment %in% c("Inulin_3", "Iron", "Control")),
    formula = formula,
    nrep = 999,
    strata = strata
  ) %>%
    bind_rows(.id = "Distance") %>%
    filter(!Distance %in% c("bray", "d_0", "d_0.5")) -> adonis_tmp
  
  m = "PCoA"
  ps %>%
    phyloseq_plot_bdiv(dlist,
                       m = m,
                       seed = 123,
                       axis1 = 1,
                       axis2 = 2) -> plots
  
  # removing some, otherwise it is too much...
  plots$bray = NULL
  plots$d_0 = NULL
  plots$d_0.5 = NULL
  
  plots %>%
    plyr::ldply(function(x) x$data) -> df
  
  names(df)[1] <- "distance"
  p = ggplot(df, aes_string(colnames(df)[2], colnames(df)[3]))
  p = p + geom_point(size=2,
                     aes_string(color= group_color, 
                                shape = group_shape,
                                alpha = group_alpha))
  
  p = p + facet_wrap( ~ distance, scales="free")
  
  p = p + ggtitle(paste0(m," using various distance metrics ")) +
    theme_light() #+ scale_color_viridis_d()
  
  p + scale_color_viridis_d() + 
    scale_alpha_continuous(range = c(0.6, 1),
                           breaks = c(20,30,40)) +
    scale_shape_manual(values = c(1,19,0,15))
  
  output = list("plot" = p,
                "PERMANOVA" = adonis_tmp)
  
  return(output)
}

###########################################################################################

betadisper <- function(dm, 
                       physeq, 
                       variable) {
  require(vegan)
  require(phyloseq)
  
  as.matrix(dm)[sample_names(physeq),sample_names(physeq)] %>%
    as.dist() -> dm
  
  
  vegan::permutest(vegan::betadisper(dm,
                                     get_variable(physeq, variable)))$tab$`Pr(>F)`[1] -> out
  
  return(out)
}

###########################################################################################

phyloseq_TW <- function(dm,
                        physeq = physeq,
                        variable = variable,
                        nrep = nrep,
                        strata = strata){
  
  as.matrix(dm)[sample_names(physeq),sample_names(physeq)] -> dm
  
  source("https://raw.githubusercontent.com/alekseyenko/WdStar/master/Wd.R")
  Tw2.posthoc.tests(dm = dm,
                    f = get_variable(physeq, variable),
                    nrep = nrep,
                    strata = get_variable(physeq, strata)) %>%
    data.frame() -> out
  out$pvalFDR = p.adjust(out$p.value, method = "fdr")
  return(out)
}
###########################################################################################

phyloseq_adonis_strata_perm <- function(dm,
                                        physeq = physeq,
                                        formula = paste0(variables, collapse=" + "),
                                        nrep = nrep,
                                        strata = strata){
  require(vegan)
  
  as.matrix(dm)[sample_names(physeq),sample_names(physeq)] %>%
    as.dist() -> dm
  
  physeq %>%
    sample_data() %>%
    data.frame() -> df
  
  if (strata %in% colnames(df)){
    perm <- how(nperm = nrep)
    setBlocks(perm) <- with(df, strata)
    
    adonis(formula = as.formula(paste("dm", paste(formula), sep=" ~ ")),
           # strata = strata,
           permutations = perm,
           data = df)$aov.tab %>% 
      data.frame() %>% 
      rownames_to_column('terms') -> out
    
    
    
  }else{
    adonis(formula = as.formula(paste("dm", paste(formula), sep=" ~ ")),
           permutations = nrep,
           data = df)$aov.tab %>% 
      data.frame() %>% 
      rownames_to_column('terms') -> out
  }
  
  
  
  return(out)
}
###########################################################################################

phyloseq_adonis <- function(dm,
                            physeq = physeq,
                            formula = paste0(variables, collapse=" + "),
                            nrep = nrep,
                            strata = strata){
  require(vegan)
  
  as.matrix(dm)[sample_names(physeq),sample_names(physeq)] %>%
    as.dist() -> dm
  
  physeq %>%
    sample_data() %>%
    data.frame() -> df
  
  if (strata %in% colnames(df)){
    
    adonis(formula = as.formula(paste("dm", paste(formula), sep=" ~ ")),
           strata = strata,
           permutations = nrep,
           data = df)$aov.tab %>% 
      data.frame() %>% 
      rownames_to_column('terms') -> out
    
    
    
  }else{
    adonis(formula = as.formula(paste("dm", paste(formula), sep=" ~ ")),
           permutations = nrep,
           data = df)$aov.tab %>% 
      data.frame() %>% 
      rownames_to_column('terms') -> out
  }
  
  
  
  return(out)
}

###########################################################################################

phyloseq_run_compare_means <- function(tmp = tmp,
                                       group = group,
                                       comp = comp,
                                       prev = prev,
                                       varcoef = varcoef)
{
  out=NULL
  
  out <- vector("list", length(tmp %>% 
                                 get_variable(group) %>% levels()))
  names(out) <- tmp %>% 
    get_variable(group) %>% levels()
  
  for(tp in tmp %>% 
      get_variable(group) %>%
      unique()){
    # print(tp)
    prune_samples(get_variable(tmp, group) == tp,
                  tmp) %>%
      transform_sample_counts(function(x) x/sum(x) * 100) %>%
      filter_taxa(function(x) sum(x > 0) > (prev*length(x)), TRUE) %>%
      filter_taxa(function(x) sd(x)/mean(x) > varcoef, TRUE) %>%
      psmelt() -> tmp2
    
    ggpubr::compare_means(formula = as.formula(paste0("Abundance ~ ", paste0(comp))),
                          group.by = c("Species"),
                          data = tmp2,
                          method = "wilcox.test",
                          p.adjust.method = "fdr") %>%
      select(-.y., -p.format, -p.signif) %>%
      arrange(p.adj) %>%
      mutate(signif = ifelse(p.adj <= 0.05, 'SIGN', 'NS')) %>% 
      mutate(Site = tp) -> results
    
    out[[tp]] <- results
  }
  return(out)
}

###########################################################################################

phyloseq_run_Deseq <- function(tmp = tmp,
                               group = group,
                               comp = comp,
                               prev = prev,
                               varcoef = varcoef,
                               A = A,
                               vsB = vsB)
{
  out=NULL
  out <- vector("list", length(tmp %>% 
                                 get_variable(group) %>% 
                                 levels()))
  names(out) <- tmp %>% 
    get_variable(group) %>% levels()
  
  for(tp in tmp %>% 
      get_variable(group) %>%
      unique()){
    # print(tp)
    prune_samples(get_variable(tmp, group) == tp,
                  tmp) %>%
      # transform_sample_counts(function(x) x/sum(x) * 100) %>%
      filter_taxa(function(x) sum(x > 0) > (prev*length(x)), TRUE) %>%
      filter_taxa(function(x) sd(x)/mean(x) > varcoef, TRUE) %>%
      phyloseq_to_deseq2(as.formula(paste0("~  ", paste0(comp)))) -> cds # convert 
    
    # calculate geometric means prior to estimate size factors
    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    geoMeans = apply(counts(cds), 1, gm_mean)
    cds = estimateSizeFactors(cds, geoMeans = geoMeans)
    
    # run DESeq function
    dds <- DESeq(cds,
                 fitType = "local") #local
    
    
    
    results <- results(dds,
                       contrast = c(comp, A, vsB),
                       tidy = TRUE) %>% 
      dplyr::rename('OTU' = row) %>%
      dplyr::rename(p.adj = padj) %>%
      arrange(p.adj) %>%
      mutate(signif = ifelse(p.adj <= 0.05, 'SIGN', 'NS')) %>% 
      mutate(Site = tp)
    
    out[[tp]] <- results
    
  }
  
  return(out)
}

###########################################################################################

phyloseq_run_ALDEx2 <- function(tmp = tmp,
                                group = group,
                                comp = comp,
                                prev = prev,
                                varcoef = varcoef,
                                mc = mc,
                                denom = denom)
{
  require(ALDEx2)
  out=NULL
  out <- vector("list", length(tmp %>% 
                                 get_variable(group) %>% 
                                 levels()))
  names(out) <- tmp %>% 
    get_variable(group) %>% levels()
  
  print(out)
  
  for(tp in tmp %>% 
      get_variable(group) %>%
      unique()){
    
    print(tp)
    
    prune_samples(get_variable(tmp, group) == tp,
                  tmp) %>%
      # transform_sample_counts(function(x) x/sum(x) * 100) %>%
      filter_taxa(function(x) sum(x > 0) > (prev*length(x)), TRUE) %>%
      filter_taxa(function(x) sd(x)/mean(x) > varcoef, TRUE) -> cds
    
    
    ALDEx2::aldex.clr(data.frame(phyloseq::otu_table(cds)), phyloseq::sample_data(cds) %>% data.frame() %>% pull(comp), 
                      mc.samples = mc, 
                      denom = denom,
                      verbose = F, useMC = TRUE) -> x
    
    x %>%
      ALDEx2::aldex.kw(useMC = TRUE) -> aldex2_da
    
    aldex2_da %>%
      rownames_to_column(var = "Species") %>%
      #filter(glm.eBH < 0.05) %>%
      #arrange(glm.eBH) %>%
      dplyr::mutate("site" = tp) %>%
      left_join(tax_table(tmp) %>% data.frame() %>%    
                  rownames_to_column('OTU'),
                by = "Species",
                suffix = c("", ".y")) -> aldex2_tax
    
    out[[tp]] <- aldex2_tax
    
  }
  return(out)
  
}

###########################################################################################

