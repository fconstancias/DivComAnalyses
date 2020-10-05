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
  detach("package:ape", unload=TRUE); detach("package:abdiv", unload=TRUE); detach("package:GUniFrac", unload=TRUE)

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


