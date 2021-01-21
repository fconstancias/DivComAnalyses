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

phyloseq_compute_bdiv <- function(phylo_tmp,
                                  norm = "pc",
                                  phylo = FALSE,
                                  seed = 123)
{
  require(tidyverse)
  require(ape)
  require(phyloseq)
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
  detach("package:ape", unload=TRUE); detach("package:GUniFrac", unload=TRUE)

}

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
#'

phyloseq_plot_bdiv <- function(ps_rare,
                               dlist,
                               m = "PCoA",
                               seed = 123,
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
          p <- p + ggtitle(paste0(m," using distance method ",   i, "\n",
                                  " NMDS 2d stress = ", iMDS$grstress %>% round(2))) +
            geom_point(size = 4) + theme_bw() #+
          # ggrepel::geom_text_repel(cex=2.5,aes(label=sample))
          plot_list[[i]] = p
        }
        if(m == "PCoA")
        {
          p <- p + ggtitle(paste0(m," using distance method ",  i)) +
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
#'
phyloseq_plot_PCoA_3d <- function(ps_rare,
                                  dlist,
                                  m = "PCoA",
                                  seed = 123,
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

calc_pairwise_permanovas_strata <- function(dm, metadata_map, compare_header, n_perm, strat) {
  # require(mctoolsr)

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
    dm_w_map_filt = filter_dm(dm_w_map, filter_cat = "in_pair",
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

  # detach("package:mctoolsr", unload=TRUE)
  return(results)
}

filter_dm <- function (input_dm, filter_cat, filter_vals, keep_vals)
{
  map_filt = test_filt_map(input_dm$map_loaded, filter_cat, filter_vals,
                       keep_vals)
  dm = as.matrix(input_dm$dm_loaded)
  samplesToUse = intersect(colnames(dm), row.names(map_filt))
  dm_use = as.dist(dm[match(samplesToUse, colnames(dm)), match(samplesToUse,
                                                               colnames(dm))])
  map_use = map_filt[match(samplesToUse, row.names(map_filt)),
  ]
  list(dm_loaded = dm_use, map_loaded = map_use)
}

test_filt_map = function(map, filter_cat, filter_vals, keep_vals){
  if(!missing(filter_vals) & !missing(keep_vals)){
    stop('Can only handle filter_vals or keep_vals, not both.')
  }
  if(!filter_cat %in% names(map)){
    stop('filter_cat not found in mapping file headers. Check spelling.')
  }
  # filter out values within mapping category
  else if(!missing(filter_cat) & !missing(filter_vals)){
    map_f = map[!map[, filter_cat] %in% filter_vals, , drop = FALSE]
    map_f = droplevels(map_f)
    if(nrow(map_f) == 0){
      stop('All rows filtered out. Check spelling of filter parameters.')
    }
  }
  # keep certain values with mapping category
  else if(!missing(filter_cat) & !missing(keep_vals)){
    map_f = map[map[,filter_cat] %in% keep_vals, , drop = FALSE]
    map_f = droplevels(map_f)
    if(nrow(map_f) == 0){
      stop('All rows filtered out. Check spelling of filter parameters.')
    }
  }
  map_f
}

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


physeq_pairwise_permanovas <- function(dm, physeq, compare_header, n_perm, strat) {
 # require(mctoolsr)

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
    dm_w_map_filt = filter_dm(dm_w_map, filter_cat = "in_pair",
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

  #detach("package:mctoolsr", unload=TRUE)

  return(results)
}

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

physeq_betadisper <- function(dm,
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

phyloseq_TW <- function(dm,
                        physeq = physeq,
                        variable = variable,
                        nrep = nrep,
                        strata = strata){

  as.matrix(dm)[sample_names(physeq),sample_names(physeq)] -> dm

  # source("https://raw.githubusercontent.com/alekseyenko/WdStar/master/Wd.R")
  Tw2.posthoc.tests(dm = dm,
                    f = get_variable(physeq, variable),
                    nrep = nrep,
                    strata = get_variable(physeq, strata)) %>%
    data.frame() -> out
  out$pvalFDR = p.adjust(out$p.value, method = "fdr")
  return(out)
}

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

    vegan::adonis(formula = as.formula(paste("dm", paste(formula), sep=" ~ ")),
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
#'plot_list %>%
#'  phyloseq_plot_ordinations_facet(color_group = "SampleType",
#'                                 shape_group = NULL)
#'


phyloseq_plot_ordinations_facet <- function(plot_list,
                                            color_group,
                                            shape_group = NULL)
{
  plot_list %>%
  plyr::ldply(function(x) x$data) -> df

names(df)[1] <- "distance"

df %>%
  ggplot(aes_string(colnames(df)[2], colnames(df)[3])) -> p

p = p + geom_point(size=2,
                   aes_string(color= color_group, 
                              shape = shape_group))

p = p + facet_wrap( ~ distance, scales="free")

p = p + ggtitle(paste0("Ordination using various distance metrics ")) +
  theme_light() 

return(p)
}

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
#'

phyloseq_ordinations_expl_var <- function(plot_list)
{
  plot_list %>%
    plyr::ldply(function(x) x$labels$x) %>%
    bind_cols(plyr::ldply(plot_list, function(x) x$labels$y)) %>%
    data.frame() -> df
  
  return(df)
}

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
#'

phyloseq_distance_boxplot <- function(p = ps, dist = dlist$wjaccard, d = "SampleType") 
  {
  
  require("phyloseq")
  require("tidyverse")
  
  s <- sample_names(p)
  as.matrix(dist)[s,s] %>%
    as.dist() -> dm
  
  # calc distances
  # wu = phyloseq::distance(p, m)
  wu.m = reshape2::melt(as.matrix(dist))
  
  # remove self-comparisons
  wu.m = wu.m %>%
    filter(as.character(Var1) != as.character(Var2)) %>%
    mutate_if(is.factor, as.character)
  
  # get sample data (S4 error OK and expected)
  sd = sample_data(p) %>%
    data.frame() %>%
    rownames_to_column("tmp") %>%
    select(tmp, all_of(d)) %>%
    mutate_if(is.factor,as.character)
  
  # combined distances with sample data
  colnames(sd) = c("Var1", "Type1")
  wu.sd = left_join(wu.m, sd, by = "Var1")
  
  colnames(sd) = c("Var2", "Type2")
  wu.sd = left_join(wu.sd, sd, by = "Var2")
  
  # plot
  p = ggplot(wu.sd, aes(x = Type2, y = value)) +
    theme_bw() +
    geom_boxplot(aes(color = ifelse(Type1 == Type2, "red", "black"))) +
    geom_jitter(aes(color = ifelse(Type1 == Type2, "red", "black")))+
    scale_color_identity() +
    facet_wrap(~ Type1, scales = "free_x") +
    theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    # ggtitle(paste0("Distance Metric = ")) +
    ylab("Distance")# +
  # xlab()
  
  # return
  out = list(plot=p,
             matrix = wu.sd)
  return(out)
  
}



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
              
phyloseq_dbRDA <- function(ps,
                           dist,
                           forumla = paste0(variables, collapse=" + "))
{
require(plyr); require(ggvegan)

  
ps %>% sample_data() %>% data.frame() -> metadata

dbRDA <- vegan::capscale(formula(paste0("dist","~",forumla)), 
                         metadata,
                         add = TRUE)

# overll significance of the model
anova(dbRDA) %>%
  data.frame() -> anova_all

# significance of different covariables
anova(dbRDA, by = "terms") %>%
  data.frame() -> anova_terms

# source('https://raw.githubusercontent.com/fawda123/ggord/master/R/ggord.R')

# ggord(dbRDA, grp_in = metadata[,variables]) -> p
autoplot(dbRDA) -> p

return(out <- list("plot" = p,
                   "dbRDA" = dbRDA,
                   "anova_all" = anova_all,
                   "anova_terms" = anova_terms))

detach("package:plyr", unload=TRUE);detach("package:ggvegan", unload=TRUE)

}


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

phyloseq_plot_dbrda <- function(physeq, dm, grouping_column, pvalueCutoff = 0.01, norm_method = NULL, 
                              env.variables = NULL, num.env.variables = NULL, exclude.variables = NULL, 
                              draw_species = F) 
{
  abund_table <- otu_table(physeq)
  meta_table <- data.frame(sample_data(physeq))[,c(env.variables,grouping_column)]
  complete.cases(meta_table) -> cc
  as.matrix(dm)[cc,cc] -> dm
  meta_table[cc,] -> meta_table
  abund_table.adonis <- adonis(formula = as.formula(paste("dm", paste(.), sep=" ~ ")),
                               permutations = 999,
                               data = meta_table[,c(env.variables)])
  bestEnvVariables <- rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)" <= 
                                                             pvalueCutoff]
  
  # abund_table.adonis <- vegan::adonis(dist_ps ~ ., data = meta_table)
  
  bestEnvVariables <- bestEnvVariables[!is.na(bestEnvVariables)]
  if (!is.null(env.variables) && (env.variables %in% bestEnvVariables)) {
    bestEnvVariables <- env.variables
  }
  if (!is.null(num.env.variables)) {
    if (num.env.variables > length(bestEnvVariables)) {
      stop(cat(paste("Choose a number less than", length(bestEnvVariables))))
    }
    else {
      bestEnvVariables <- bestEnvVariables[1:num.env.variables]
    }
  }
  if (!is.null(exclude.variables) && (exclude.variables %in% 
                                      bestEnvVariables)) {
    bestEnvVariables <- bestEnvVariables[!(bestEnvVariables %in% 
                                             exclude.variables)]
  }
  eval(parse(text = paste("sol <- capscale(dm ~ ", do.call(paste, 
                                                           c(as.list(bestEnvVariables), sep = " + ")), ",data=meta_table)", 
                          sep = "")))
  scrs <- vegan::scores(sol, display = c("sp", "wa", "lc", 
                                         "bp", "cn"))
  df_sites <- data.frame(scrs$sites, meta_table[, grouping_column])
  colnames(df_sites) <- c("x", "y", "Groups")
  p <- ggplot2::ggplot()
  p <- p + ggplot2::geom_point(data = df_sites, aes(x, y, colour = Groups))
  multiplier <- vegan:::ordiArrowMul(scrs$biplot)
  df_arrows <- scrs$biplot * multiplier
  colnames(df_arrows) <- c("x", "y")
  df_arrows = as.data.frame(df_arrows)
  p <- p + geom_segment(data = df_arrows, aes(x = 0, y = 0, 
                                              xend = x, yend = y), arrow = arrow(length = unit(0.2, 
                                                                                               "cm")), color = "#808080", alpha = 0.5)
  p <- p + geom_text(data = as.data.frame(df_arrows * 1.1), 
                     aes(x, y, label = rownames(df_arrows)), color = "#808080", 
                     alpha = 0.5)
  df_species <- as.data.frame(scrs$species)
  colnames(df_species) <- c("x", "y")
  if (draw_species) {
    p <- p + geom_point(data = df_species, aes(x, y, shape = "Species")) + 
      scale_shape_manual("", values = 2)
  }
  p <- p + theme_bw() + xlab("CCA1") + ylab("CCA2")
  return(list(plot= p,
              capscale = sol))
}