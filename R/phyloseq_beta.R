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
#'esophagus %>%
#'phyloseq_plot_bdiv(., dlist = dist) -> ords
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
                               axis2 = axis2,
                               TSNE_per = 5)
{
  if (m == "CoDa")
  {
    ps_rare <- microbiome::transform(ps_rare, "clr")

    ord_clr <- phyloseq::ordinate(ps_rare, "RDA")

    #Scale axes and plot ordination
    phyloseq::plot_ordination(ps_rare, ord_clr, type="samples") -> p
    # coord_fixed(ord_clr$CA$eig[2] / sum(ord_clr$CA$eig) / ord_clr$CA$eig[1] / sum(ord_clr$CA$eig))
    #https://github.com/ggloor/CoDa_microbiome_tutorial/wiki/Part-1:-Exploratory-Compositional-PCA-biplot


    aidist <- ps_rare %>%
      phyloseq::distance(method = "euclidean")

    out <- list("PCA" = p,
                "physeq_clr" = ps_rare,
                "aidist" = aidist)

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
        tsne_out <- Rtsne::Rtsne(dlist[[i]], dims = 2, perplexity = TSNE_per, verbose = T)
        proj <- tsne_out$Y %>% data.frame()

        rownames(proj) <- rownames(t(otu_table(ps_rare)))

        proj2 <- cbind(proj, sample_data(ps_rare))

        # rownames(proj) == sample_data(ps_rare)$SampleID

        # microbiome::plot_landscape(proj, legend = T, size = 1)
        p <- phyloseq::plot_ordination(ps_rare,
                                       proj)

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
          stress_list <- vector("list", length(dlist))
          names(stress_list) =  names(dlist)

          stress = iMDS$grstress %>% round(2)

          p <- p + ggtitle(paste0(m," using distance method ",   i, "\n",
                                  " NMDS 2d stress = ", stress)) +
            geom_point(size = 4) + theme_bw() #+
          # ggrepel::geom_text_repel(cex=2.5,aes(label=sample))
          plot_list[[i]] = p

          stress_list[[i]] = stress

          plot_list <- c(plot_list, stress_list)
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

calc_pairwise_permanovas_strata <- function(dm, physeq, compare_header, n_perm,  strat) {
  # require(mctoolsr)

  physeq %>%
    sample_data() %>%
    data.frame() -> metadata_map

  as.matrix(dm)[sample_names(physeq),sample_names(physeq)] %>%
    as.dist() -> dm

  comp_var = as.factor(metadata_map[, compare_header])
  comp_pairs = utils::combn(base::levels(comp_var), 2)
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
                                                            strat], na.action = na.exclude)
      }
      else {
        m = vegan::adonis(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[,
                                                                             compare_header],
                          strata = dm_w_map_filt$map_loaded[,
                                                            strat],
                          na.action = na.exclude)
      }
    }else{
      if (!missing(n_perm)) {
        m = vegan::adonis(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[,
                                                                             compare_header], permutations = n_perm,
                          na.action = na.exclude)
      }
      else {
        m = vegan::adonis(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[,
                                                                             compare_header] %>% as.factor(),
                          na.action = na.exclude)
      }
    }
    pval = c(pval, m$`Pr(>F)`[1])
    R2 = c(R2, m$R2[1])
  }
  results = data.frame(t(comp_pairs), R2, pval)
  results$pvalBon = pval * length(pval)
  results$pvalFDR = round(pval * (length(pval)/rank(pval, ties.method = "average")),
                          3)

  # detach("package:mctoolsr", unload=TRUE)
  return(results)

  detach("package:vegan", unload=TRUE)

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
#'source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_beta.R")
#'library(phyloseq); library(tidyverse)
#'data(enterotype)
#'phyloseq_compute_bdiv(enterotype) -> dist
#'enterotype %>%
#'subset_samples(!is.na(Nationality)) %>%
#'physeq_pairwise_permanovas(dm = dist$bray, physeq = ., compare_header = "Nationality")
#'or
#'physeq_pairwise_permanovas_adonis2(dm = dist$bray, physeq = ., compare_header = "Nationality")
#'


physeq_pairwise_permanovas <- function(dm, physeq, compare_header, n_perm, strat = FALSE, terms_margins = "terms") {
  require(vegan)

  as.matrix(dm)[sample_names(physeq),sample_names(physeq)] %>%
    as.dist() -> dm

  physeq %>%
    sample_data() %>%
    data.frame() -> metadata_map

  comp_var = as.factor(metadata_map[, compare_header])
  comp_pairs = utils::combn(levels(comp_var), 2)
  pval = c()
  R2 = c()
  for (i in 1:ncol(comp_pairs)) {
    pair = comp_pairs[, i]
    dm_w_map = base::list(dm_loaded = dm, map_loaded = metadata_map)
    dm_w_map$map_loaded$in_pair = comp_var %in% pair
    dm_w_map_filt = filter_dm(dm_w_map, filter_cat = "in_pair",
                              keep_vals = TRUE)

    if (strat %in% colnames(metadata_map)){

      if (!missing(n_perm)) {
        m = vegan::adonis2(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[,
                                                                              compare_header], permutations = n_perm,
                           strata = dm_w_map_filt$map_loaded[,
                                                             strat],
                           by = terms_margins)
      }
      else {
        m = vegan::adonis2(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[,
                                                                              compare_header],
                           strata = dm_w_map_filt$map_loaded[,
                                                             strat],
                           by = terms_margins)
      }
    }else{
      if (!missing(n_perm)) {
        m = vegan::adonis2(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[,
                                                                              compare_header], permutations = n_perm,
                           by = terms_margins)
      }
      else {
        m = vegan::adonis2(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[,
                                                                              compare_header],
                           by = terms_margins)
      }
    }
    pval = c(pval, m$`Pr(>F)`[1])
    R2 = c(R2, m$R2[1])
  }
  results = data.frame(t(comp_pairs), R2, pval)
  results$pvalBon = pval * length(pval)
  results$pvalFDR = round(pval * (length(pval)/rank(pval, ties.method = "average")),
                          3)

  #detach("package:mctoolsr", unload=TRUE)

  return(results)



  detach("package:vegan", unload=TRUE)

}


physeq_pairwise_permanovas_adonis2 <- function(dm, physeq, compare_header, n_perm = 999, strata = "none", terms_margins = "terms") {

  require(vegan)

  as.matrix(dm)[sample_names(physeq),sample_names(physeq)] %>%
    as.dist() -> dist

  physeq %>%
    sample_data() %>%
    data.frame() -> df

  comp_var = as.factor(df[, compare_header])
  comp_pairs = combn(levels(comp_var), 2)

  pval= NULL; R2 = NULL

  for (i in 1:ncol(comp_pairs)) {

    pair = comp_pairs[, i]
    dm_w_map = list(dm_loaded = dm, map_loaded = df)
    dm_w_map$map_loaded$in_pair = comp_var %in% pair
    dm_w_map_filt = filter_dm(dm_w_map, filter_cat = "in_pair",
                              keep_vals = TRUE)

    dm_w_map_filt$dm_loaded -> dist_tmp
    dm_w_map_filt$map_loaded -> df_tmp

    if (strata %in% colnames(df)){

      df_tmp %>%
        mutate("strata" = get(strata)) -> df_tmp

      perm <- how(nperm = n_perm)
      setBlocks(perm) <- with(df_tmp, strata)

      adonis2(formula = as.formula(paste("dist_tmp", paste(compare_header), sep=" ~ ")),
              permutations = perm,
              by = terms_margins,
              data = df_tmp) %>%
        as.data.frame() -> m

    }else{

      adonis2(formula = as.formula(paste("dist_tmp", paste(compare_header), sep=" ~ ")),
              permutations = n_perm,
              by = terms_margins,
              data = df_tmp) %>%
        as.data.frame() -> m
    }
    pval = c(pval, m$`Pr(>F)`[1])
    R2 = c(R2, m$R2[1])
  }
  results = data.frame(t(comp_pairs), R2, pval)
  results$pvalBon = pval * length(pval)
  results$pvalFDR = round(pval * (length(pval)/rank(pval, ties.method = "average")),
                          3)

  #detach("package:mctoolsr", unload=TRUE)

  return(results)


  detach("package:vegan", unload=TRUE)

}


filter_dm <- function (input_dm, filter_cat, filter_vals, keep_vals){
  map_filt = test_filt_map(input_dm$map_loaded, filter_cat, filter_vals,
                           keep_vals)
  dm = as.matrix(input_dm$dm_loaded)
  samplesToUse = base::intersect(colnames(dm), row.names(map_filt))
  dm_use = as.dist(dm[match(samplesToUse, colnames(dm)), match(samplesToUse,
                                                               colnames(dm))])
  map_use = map_filt[match(samplesToUse, row.names(map_filt)),
  ]
  base::list(dm_loaded = dm_use, map_loaded = map_use)
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
    map_f = base::droplevels(map_f)
    if(nrow(map_f) == 0){
      stop('All rows filtered out. Check spelling of filter parameters.')
    }
  }
  # keep certain values with mapping category
  else if(!missing(filter_cat) & !missing(keep_vals)){
    map_f = map[map[,filter_cat] %in% keep_vals, , drop = FALSE]
    map_f = base::droplevels(map_f)
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
                                     get_variable(physeq, variable)))$tab$`Pr(>F)`[1] -> pval

  boxplot(vegan::betadisper(dm,
                            get_variable(physeq, variable)),las=2,
          main=paste0("Multivariate Dispersion Test "," pvalue = ",
                      vegan::permutest(betadisper(bc, get_variable(physeq, variable)))$tab$`Pr(>F)`[1])) -> plot

  return(out <- list("pval" = pval,
                     "plot" = plot))


  detach("package:vegan", unload=TRUE)

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
                        nrep = 999,
                        strata = NULL){

  as.matrix(dm)[sample_names(physeq),sample_names(physeq)] -> dm

  # source("https://raw.githubusercontent.com/alekseyenko/WdStar/master/Wd.R")
  if(!is.null(strata)){strata = get_variable(physeq, strata)}

  Tw2.posthoc.tests(dm = dm,
                    f = get_variable(physeq, variable),
                    nrep = nrep,
                    strata = strata) %>%
    # as.data.frame.table() %>%
    as.data.frame.array() %>%
    mutate(pvalFDR = p.adjust(p.value, method = "fdr")) -> out

  return(out)
}

is.dist = function(x) any(class(x)=='dist')

dist.sigma2 = function(dm){
  dd = as.matrix(dm)
  dd[upper.tri(dd)]=0 ##
  sum(dd^2)/nrow(dd)/(nrow(dd)-1)
}

dist.ss2 = function(dm2, f){ #dm2 is matrix of square distances; f factor
  K = sapply(levels(f), function(lev) f==lev)
  t(K)%*%dm2%*%K/2
}

dist.group.sigma2 = function(dm, f){
  diag(dist.ss2(as.matrix(dm)^2, f))/table(f)/(table(f)-1)
}

dist.cohen.d = function(dm, f){
  if(nlevels(f) != 2) return(NULL)
  SS2 = dist.ss2(as.matrix(dm^2),f)
  ns = summary(f)
  N = sum(ns)

  SST = sum(SS2)/N
  SSW = sum(diag(SS2)/ns)

  mean.diff = (sqrt((ns[1]+ns[2])/(ns[1]*ns[2])))*sqrt(SST-SSW)

  sigmas = diag(SS2)/ns/(ns-1)
  s1 = sigmas[1]
  s2 = sigmas[2]

  mean.diff/sqrt(((ns[1]-1)*s1 + (ns[2]-1)*s2)/(sum(ns)-2))
}

Tw2 = function(dm, f){
  if(nlevels(f) != 2) return(NULL)
  SS2 = dist.ss2(as.matrix(dm^2),f)
  ns = summary(f)
  N = sum(ns)

  SST = sum(SS2)/N
  SSW1 = SS2[1,1]/ns[1]
  SSW2 = SS2[2,2]/ns[2]
  SSW = SSW1 + SSW2

  s1 = SSW1/(ns[1]-1)
  s2 = SSW2/(ns[2]-1)

  t.stat = (ns[1]+ns[2])/(ns[1]*ns[2])*(SST-SSW)/(s1/ns[1] + s2/ns[2])
  t.stat
}

WdS = function(dm, f){
  # This method computes Wd* statistic for distance matrix dm and factor f
  ns = table(f)
  SS2 = dist.ss2(as.matrix(dm)^2, f)
  s2 = diag(SS2)/ns/(ns-1)
  W = sum(ns/s2)

  idxs = apply(combn(levels(f), 2),2, function(idx) levels(f) %in% idx)

  Ws = sum(apply(idxs, 2,
                 function(idx) sum(ns[idx])/prod(s2[idx]) *
                   (sum(SS2[idx, idx])/sum(ns[idx]) - sum(diag(SS2[idx, idx])/ns[idx]))))
  k=nlevels(f)
  h = sum( (1-ns/s2/W)^2/(ns-1))
  Ws/W/(k-1)/(1+(2*(k-2)/(k^2-1))*h)
}

generic.distance.permutation.test =
  function(test.statistic, dm, f, nrep=999, strata = NULL){
    N = length(f)
    generate.permutation=function(){
      f[sample(N)]
    }

    if(!is.null(strata)){
      # map elements of each strata back to their positions in the factor variable
      strata.map = order(unlist(tapply(seq_along(f), strata, identity)))
      generate.permutation=function(){
        p = unlist(tapply(f,strata,sample)) # permute within strata
        p[strata.map]
      }
    }

    stats = c(test.statistic(dm, f),
              replicate(nrep,
                        test.statistic(dm, generate.permutation())))

    p.value = sum(stats>=stats[1])/(nrep+1)
    statistic = stats[1]
    list(p.value = p.value, statistic = statistic, nrep=nrep)
  }

Tw2.test = function(dm, f, nrep=999, strata=NULL){
  generic.distance.permutation.test(Tw2, dm = dm, f = f, nrep = nrep, strata=strata)
}

WdS.test = function(dm, f, nrep=999, strata=NULL){
  generic.distance.permutation.test(WdS, dm = dm, f = f, nrep = nrep, strata=strata)
}

Tw2.posthoc.tests = function(dm, f, nrep=999, strata=NULL){
  dd = as.matrix(dm)
  Tw2.subset.test=function(include.levels){
    subs = f %in% include.levels
    c(include.levels,
      table(f[subs, drop=T]),
      Tw2.test(dd[subs, subs], f[subs,drop=T], nrep=nrep, strata=strata[subs]))
  }
  res = t(combn(levels(f), 2, Tw2.subset.test))
  colnames(res) = c("Level1", "Level2", "N1", "N2", "p.value", "tw2.stat", "nrep")
  res
}

Tw2.posthoc.1vsAll.tests = function(dm, f, nrep=999, strata=NULL){
  Tw2.subset.test=function(level){
    fs = factor(f == level)
    c(table(fs), Tw2.test(dm, fs, nrep=nrep, strata=strata))
  }
  res = t(sapply(levels(f), Tw2.subset.test))
  colnames(res) = c("N1", "N2", "p.value", "tw2.stat", "nrep")
  res
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
                                        physeq,
                                        formula = paste0(variables, collapse=" + "),
                                        nrep,
                                        terms_margins = "terms",
                                        strata = "none"){
  require(vegan)

  as.matrix(dm)[sample_names(physeq),sample_names(physeq)] %>%
    as.dist() -> dm

  physeq %>%
    sample_data() %>%
    data.frame() %>%
    rownames_to_column("tmp_id") -> df

  if (strata %in% colnames(df)){

    df %>%
      mutate("strata" = get(strata)) -> df

    perm <- how(nperm = nrep)
    setBlocks(perm) <- with(df, strata)

    adonis2(formula = as.formula(paste("dm", paste(formula), sep=" ~ ")),
            # strata = strata,
            permutations = perm,
            by = terms_margins,
            data = df) %>%
      data.frame() %>%
      rownames_to_column('terms') -> out



  }else{
    adonis2(formula = as.formula(paste("dm", paste(formula), sep=" ~ ")),
            permutations = nrep,
            by = terms_margins,
            data = df) %>%
      data.frame() %>%
      rownames_to_column('terms') -> out
  }



  return(out)

  detach("package:vegan", unload=TRUE)

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
#'data(enterotype)
#'phyloseq_compute_bdiv(enterotype) -> dist
#'enterotype %>%
#'subset_samples(!is.na(Gender)) %>%
#'phyloseq_adonis2(dm = dist$bray, physeq = ., formula = "Gender")


phyloseq_adonis2 <- function(dm,
                             physeq = physeq,
                             formula = paste0(variables, collapse=" + "),
                             nrep = 999,
                             strata = "none",
                             terms_margins = "terms"){
  require(vegan)

  as.matrix(dm)[sample_names(physeq),sample_names(physeq)] %>%
    as.dist() -> dm

  physeq %>%
    sample_data() %>%
    data.frame() -> df

  if (strata %in% colnames(df)){

    adonis2(formula = as.formula(paste("dm", paste(formula), sep=" ~ ")),
            strata = strata,
            permutations = nrep,
            data = df,
            by = terms_margins) %>%
      data.frame() %>%
      rownames_to_column('terms') -> out



  }else{
    adonis2(formula = as.formula(paste("dm", paste(formula), sep=" ~ ")),
            permutations = nrep,
            data = df,
            by = terms_margins) %>%
      data.frame() %>%
      rownames_to_column('terms') %>%
      dplyr::rename('Pr(>F)' = `Pr..F.` ) -> out
  }



  return(out)

  detach("package:vegan", unload=TRUE)

}

phyloseq_adonis <- function(dm,
                            physeq = physeq,
                            formula = paste0(variables, collapse=" + "),
                            nrep = 999,
                            strata = "none",
                            terms_margins = "terms"){
  require(vegan)

  as.matrix(dm)[sample_names(physeq),sample_names(physeq)] %>%
    as.dist() -> dm

  physeq %>%
    sample_data() %>%
    data.frame() -> df

  if (strata %in% colnames(df)){

    adonis2(formula = as.formula(paste("dm", paste(formula), sep=" ~ ")),
            strata = strata,
            permutations = nrep,
            data = df,
            by = terms_margins) %>%
      data.frame() %>%
      rownames_to_column('terms') -> out



  }else{
    adonis2(formula = as.formula(paste("dm", paste(formula), sep=" ~ ")),
            permutations = nrep,
            data = df,
            by = terms_margins) %>%
      data.frame() %>%
      rownames_to_column('terms') -> out
  }



  return(out)

  detach("package:vegan", unload=TRUE)

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
                                            color_group = "treatment_grouped",
                                            alpha = 0.9,
                                            shape_group = NULL,
                                            axis_names = c("Axis.1", "Axis.2"),
                                            return_eig_df = FALSE)
{
  ## ------------------------------------------------------------------------

  if(!any(class(plot_list[[1]]) == "ggplot")){

    ## ------------------------------------------------------------------------


    unlist(plot_list, recursive = FALSE, use.names = TRUE) -> un_listed

    un_listed %>%
      plyr::ldply(function(x) x$data) %>%
      dplyr::rename(group = `.id`) %>%
      separate(group, into = c("a","b"), sep = "\\.") ->  df

    df %>%
      ggplot(aes_string(axis_names[[1]], axis_names[[2]])) -> p
    # tidyr::separate(group, into = c("gpA", "gpB"), sep = ".", fill= "left", remove = FALSE)
    p = p + geom_point(size=2,
                       aes_string(color= color_group,
                                  shape = shape_group,
                                  alpha = alpha))

    p = p + facet_wrap(b ~ a, scales="free")

    p = p + ggtitle(paste0("Ordination using various distance metrics ")) +
      theme_light()

    ## ------------------------------------------------------------------------
    out <- p

    if(return_eig_df!=FALSE){
      un_listed %>%
        plyr::ldply(function(x) x$labels$x) %>%
        bind_cols(plyr::ldply(un_listed, function(x) x$labels$y)) %>%
        data.frame() -> eig_df

      out <- list("plots" = p,
                  "eig_df" = eig_df)
    }
  }else{

    ## ------------------------------------------------------------------------

    plot_list %>%
      plyr::ldply(function(x) x$data) %>%
      dplyr::rename(distance = `.id`) -> df

    df %>%
      ggplot(aes_string(axis_names[[1]], axis_names[[2]])) -> p
    # tidyr::separate(group, into = c("gpA", "gpB"), sep = ".", fill= "left", remove = FALSE)
    p = p + geom_point(size=2,
                       aes_string(color= color_group,
                                  shape = shape_group,
                                  alpha = alpha))

    p = p + facet_wrap(distance ~. , scales="free")

    p = p + ggtitle(paste0("Ordination using various distance metrics ")) +
      theme_light()

    out <- p
    ## ------------------------------------------------------------------------
    if(return_eig_df!=FALSE){

      plot_list %>%
        plyr::ldply(function(x) x$labels$x) %>%
        bind_cols(plyr::ldply(plot_list, function(x) x$labels$y)) %>%
        data.frame() -> eig_df

      out <- list("plots" = p,
                  "eig_df" = eig_df)
    }
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

phyloseq_distance_boxplot <- function(p, dist = dlist$wjaccard, d = "SampleType", filter = NULL)
{

  require("phyloseq")
  require("tidyverse")

  as.matrix(dist)[sample_names(p),sample_names(p)] %>%
    as.dist() -> dist

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

  if(!is.null(filter)){
    wu.sd %>%
      dplyr::filter(Type1 == !!filter) -> wu.sd
  }

  # plot
  p = ggplot(wu.sd, aes(x = Type2, y = value)) +
    theme_bw() +
    geom_jitter(aes(color = ifelse(Type1 == Type2, "red", "black")),
                alpha = 0.1)+
    geom_boxplot(aes(color = ifelse(Type1 == Type2, "red", "black")),
                 fill = NA,
                 outlier.shape = NA,
                 outlier.colour = NA) +
    # outlier.shape = NA,) +
    scale_color_identity() +
    facet_wrap(~ Type1, scales = "free_x") +
    theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    # ggtitle(paste0("Distance Metric = ")) +
    ylab("Distance") + xlab(NULL)

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

phyloseq_add_taxa_vector_fix <- function(dist =NULL,
                                         phyloseq = NULL,
                                         figure_ord = pca_treat,
                                         m = "PCoA",
                                         pval_cutoff = 0.05,
                                         top_r = 12,
                                         taxrank_glom = NULL,
                                         tax_rank_plot = taxrank_glom,
                                         vector_color = "grey10",
                                         taxnames_rm = c("unknown", "Unclassified"),
                                         fact = 3,
                                         seed = 123,
                                         perm = 999,
                                         adj_method = "fdr",
                                         ggrepel_force = 25)
{
  ####----------------------

  require(phyloseq); require(tidyverse); require(vegan)

  ####----------------------

  as.matrix(dist)[sample_names(phyloseq),sample_names(phyloseq)] %>%
    as.dist() -> dist

  ####---------------------- Calculate ordination

  set.seed(seed)
  iMDS  <- ordinate(phyloseq,
                    m,
                    dist)

  ####----------------------Normalize features

  phyloseq %>%
    transform_sample_counts(function(x) {x/sum(x)} * 100)  -> tmp1
  ####----------------------

  if(taxrank_glom != FALSE) {
    tmp1 %>%
      speedyseq::tax_glom(taxrank_glom) -> ps_glom
  }

  ps_glom %>%
    otu_table() %>%
    t() %>%
    data.frame() -> tmp

  ####----------------------Create plot, store as temp variable, p
  set.seed(seed)
  p <- phyloseq::plot_ordination(phyloseq, iMDS)

  dune.spp.fit <- envfit(iMDS$vectors, tmp, permutations = perm) # this fits species vectors

  spp.scrs <- as.data.frame(scores(dune.spp.fit, display = "vectors")) #save species intrinsic values into dataframe

  spp.scrs <- cbind(spp.scrs) #add species names to dataframe
  spp.scrs <- cbind(spp.scrs, pval = dune.spp.fit$vectors$pvals, r = dune.spp.fit$vectors$r) #add pvalues to dataframe so you can select species which are significant
  #spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
  # sig.spp.scrs <- filter(spp.scrs, pval<=pval_cutoff ) %>% top_n(top_r, r) #subset data to show species significant at 0.05
  sig.spp.scrs <- spp.scrs

  ####----------------------

  as(tax_table(ps_glom), "matrix") %>%
    as.data.frame() -> tax_table

  # if(join_cbind == "join"){
  #   left_join(sig.spp.scrs,
  #             tax_table,
  #             by = c("id" = id_taxa)) %>%
  #     dplyr::rename(tax_rank_plot = all_of(tax_rank_plot)) %>%
  #     dplyr::filter(!tax_rank_plot %in% taxnames_rm,
  #                   pval<=pval_cutoff) %>%
  #     top_n(top_r, r) -> all
  # }
  # if(join_cbind == "cbind"){
  cbind(sig.spp.scrs, tax_table) %>%
    rstatix::adjust_pvalue(p.col = "pval",
                           method = adj_method ) %>%
    dplyr::rename(tax_rank_plot = all_of(tax_rank_plot)) -> envfit_all

  envfit_all %>%
    dplyr::filter(!tax_rank_plot %in% taxnames_rm,
                  pval.adj<=pval_cutoff) %>%
    top_n(top_r, r)  -> all
  # }

  # !!variable := name_of_col_from_df

  ####----------------------

  figure_ord  +
    geom_segment(data = all,
                 aes(x = 0, xend=Axis.1* fact, y=0, yend=Axis.2 * fact), arrow = arrow(length = unit(0.25, "cm")), colour = vector_color, lwd=0.3, inherit.aes = FALSE) + #add vector arrows of significant species
    ggrepel::geom_text_repel(data = all, aes(x= Axis.1* fact, y=Axis.2*fact, label = tax_rank_plot), cex = 3, direction = "both", segment.size = 0.25, inherit.aes = FALSE, force = ggrepel_force, segment.linetype = 2, segment.color = "gray70") -> p2

  ####----------------------


  ggplot() + theme_void() +
    geom_segment(data = all,
                 aes(x = 0, xend=Axis.1* fact, y=0, yend=Axis.2 * fact), arrow = arrow(length = unit(0.25, "cm")), colour = vector_color, lwd=0.3, inherit.aes = FALSE) + #add vector arrows of significant species
    ggrepel::geom_text_repel(data = all, aes(x= Axis.1* fact, y=Axis.2*fact, label = tax_rank_plot), cex = 3, direction = "both", segment.size = 0.25, inherit.aes = FALSE, force = ggrepel_force, segment.linetype = 2, segment.color = "gray70") -> p3


  ####----------------------

  out <- list("plot" = p2,
              "vectors" = p3,
              "ord" = iMDS,
              "envfit" = envfit_all,
              "signenvfit" = all)

  return(out)

  detach("package:vegan", unload=TRUE)
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
#'
phyloseq_add_taxa_vector <- function(dist,
                                     phyloseq,
                                     figure_ord = figure_pca,
                                     m = "PCoA",
                                     pval_cutoff = 0.05,
                                     top_r = 12,
                                     taxrank_glom = "Family",
                                     tax_rank_plot = "Family",
                                     id_taxa = "ASV",
                                     taxnames_rm = "unknown",
                                     fact = 3,
                                     seed = 123,
                                     perm = 999,
                                     join_cbind = "join") # to avoid issue with rownames when specific character
{
  require(phyloseq); require(tidyverse); require(vegan)

  as.matrix(dist)[sample_names(phyloseq),sample_names(phyloseq)] %>%
    as.dist() -> dist

  # Calculate ordination
  set.seed(seed)
  iMDS  <- ordinate(phyloseq,
                    m,
                    dist)

  phyloseq %>%
    transform_sample_counts(function(x) {x/sum(x)} * 100)  -> tmp1

  #
  # if (taxrank_glom != "Strain"){
  #   tmp1 %>%
  #   speedyseq::tax_glom(taxrank = taxrank_glom) -> tmp1
  #
  #   prune_taxa(data.frame(tax_table(tmp1)[,taxrank_glom])  %>%
  #                dplyr::filter(!get(taxrank_glom) %in% taxnames_rm) %>% rownames(),tmp1) -> tmp1
  #
  #   taxa_names(tmp1) <-  tax_table(tmp1)[,taxrank_glom]
  # }

  if(taxrank_glom != FALSE) {
    tmp1 %>%
      speedyseq::tax_glom(taxrank_glom) %>%
      otu_table() %>%
      t() %>%
      data.frame() -> tmp
  }else{
    tmp1 %>%
      otu_table() %>%
      t() %>%
      data.frame() -> tmp
  }


  # Create plot, store as temp variable, p
  set.seed(seed)
  p <- phyloseq::plot_ordination(phyloseq, iMDS)

  dune.spp.fit <- envfit(iMDS$vectors, tmp, permutations = perm) # this fits species vectors

  spp.scrs <- as.data.frame(scores(dune.spp.fit, display = "vectors")) #save species intrinsic values into dataframe
  spp.scrs <- cbind(spp.scrs, id = rownames(spp.scrs)) #add species names to dataframe
  spp.scrs <- cbind(spp.scrs, pval = dune.spp.fit$vectors$pvals, r = dune.spp.fit$vectors$r) #add pvalues to dataframe so you can select species which are significant
  #spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
  # sig.spp.scrs <- filter(spp.scrs, pval<=pval_cutoff ) %>% top_n(top_r, r) #subset data to show species significant at 0.05
  sig.spp.scrs <- spp.scrs

  # left_join(sig.spp.scrs,
  #          tmp1 %>%
  #   tax_table() %>%
  #   as.data.frame(),
  #   by = c("id" = id_taxa)) -> all

  as(tax_table(tmp1), "matrix") %>%
    as.data.frame() %>%
    rownames_to_column('ASV') -> tax_table

  if(join_cbind == "join"){
    left_join(sig.spp.scrs,
              tax_table,
              by = c("id" = id_taxa)) %>%
      dplyr::rename(tax_rank_plot = all_of(tax_rank_plot)) %>%
      dplyr::filter(!tax_rank_plot %in% taxnames_rm,
                    pval<=pval_cutoff) %>%
      top_n(top_r, r) -> all
  }
  if(join_cbind == "cbind"){
    cbind(sig.spp.scrs, tax_table) %>%
      dplyr::rename(tax_rank_plot = all_of(tax_rank_plot)) %>%
      dplyr::filter(!tax_rank_plot %in% taxnames_rm,
                    pval<=pval_cutoff) %>%
      top_n(top_r, r) -> all
  }


  # !!variable := name_of_col_from_df

  figure_ord  +
    geom_segment(data = all,
                 aes(x = 0, xend=Axis.1* fact, y=0, yend=Axis.2 * fact), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3, inherit.aes = FALSE) + #add vector arrows of significant species
    ggrepel::geom_text_repel(data = all, aes(x= Axis.1* fact, y=Axis.2*fact, label = tax_rank_plot), cex = 3, direction = "both", segment.size = 0.25, inherit.aes = FALSE) -> p2

  out <- list("plot" = p2,
              "ord" = iMDS,
              "envfit" = spp.scrs,
              "signenvfit" = all)

  return(out)

  detach("package:vegan", unload=TRUE)

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
#'dist=d_list$wjaccard
#'phyloseq=physeq_rare_cec_mIMT2
#'figure_ord = pca_tmp


phyloseq_add_metadata_vector <- function(dist,
                                         phyloseq,
                                         figure_ord = figure_pca,
                                         m = "PCoA",
                                         pval_cutoff = 0.05,
                                         top_r = 12,
                                         metadata_sel = c("B24_Acetate", "B24_Propionate", "B24_Butyrate", "B24_Formate", "B24_Succinate", "B24_Lactate", "B24_BCFA"),
                                         fact = 0.5,
                                         seed = 123,
                                         perm = 999,
                                         norm_method = "center_scale",
                                         color = "green",
                                         linetype = "dashed",
                                         na.rm = TRUE)
{
  require(phyloseq); require(tidyverse); require(vegan)

  as.matrix(dist)[sample_names(phyloseq),sample_names(phyloseq)] %>%
    as.dist() -> dist

  # Calculate ordination
  set.seed(seed)

  iMDS  <- ordinate(phyloseq,
                    m,
                    dist)

  phyloseq %>%
    transform_sample_counts(function(x) {x/sum(x)} * 100)  -> tmp1

  #
  # if (taxrank_glom != "Strain"){
  #   tmp1 %>%
  #   speedyseq::tax_glom(taxrank = taxrank_glom) -> tmp1
  #
  #   prune_taxa(data.frame(tax_table(tmp1)[,taxrank_glom])  %>%
  #                dplyr::filter(!get(taxrank_glom) %in% taxnames_rm) %>% rownames(),tmp1) -> tmp1
  #
  #   taxa_names(tmp1) <-  tax_table(tmp1)[,taxrank_glom]
  # }

  tmp1 %>%
    sample_data()%>%
    data.frame() %>%
    dplyr::select(any_of(metadata_sel)) %>%
    drop_na() -> tmp


  if(nrow(tmp) == 0) {
    print ("No metadata were selected - or only NA.")}else{



      if (norm_method == "center_scale")
      {
        tmp %>%
          mutate_if(is.numeric, scale) -> tmp
      }

      # Create plot, store as temp variable, p
      set.seed(seed)
      p <- phyloseq::plot_ordination(phyloseq, iMDS)


      dune.spp.fit <- envfit(iMDS$vectors, tmp, permutations = perm, na.rm = na.rm) # this fits species vectors


      spp.scrs <- as.data.frame(scores(dune.spp.fit, display = "vectors")) #save species intrinsic values into dataframe
      spp.scrs <- cbind(spp.scrs, id = rownames(spp.scrs)) #add species names to dataframe
      spp.scrs <- cbind(spp.scrs, pval = dune.spp.fit$vectors$pvals, r = dune.spp.fit$vectors$r) #add pvalues to dataframe so you can select species which are significant
      #spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
      # sig.spp.scrs <- filter(spp.scrs, pval<=pval_cutoff ) %>% top_n(top_r, r) #subset data to show species significant at 0.05
      sig.spp.scrs <- spp.scrs


      sig.spp.scrs %>%
        rownames_to_column('label') %>%
        dplyr::filter(pval<=pval_cutoff) %>%
        top_n(top_r, r) -> all



      # !!variable := name_of_col_from_df

      figure_ord  +
        geom_segment(data = all,
                     aes(x = 0, xend=Axis.1* fact, y=0, yend=Axis.2 * fact), arrow = arrow(length = unit(0.25, "cm")),linetype= linetype, colour = color, lwd=0.3, inherit.aes = FALSE) + #add vector arrows of significant species
        ggrepel::geom_text_repel(data = all, aes(x= Axis.1* fact, y=Axis.2*fact, label = label), cex = 3, direction = "both", colour = color, segment.size = 0.25, inherit.aes = FALSE) -> p2

      out <- list("plot" = p2,
                  "ord" = iMDS,
                  "envfit" = spp.scrs,
                  "signenvfit" = all)

      return(out)
    }

  detach("package:vegan", unload=TRUE)

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
                           dm,
                           forumla = paste0(variables, collapse=" + "),
                           group_plot,
                           vec_ext = 0.2)

  #TODO: https://rstudio-pubs-static.s3.amazonaws.com/694016_e2d53d65858d4a1985616fa3855d237f.html
  #   require(ggordiplots)
  # #  devtools::install_github("jfq3/ggordiplots")
  #
  # gg_ordiplot(dbRDA, groups = metadata$Group, pt.size = 3)


{
  ### ------
  require(ggvegan); require(ggord); require(vegan);

  as.matrix(dm)[sample_names(ps),sample_names(ps)] %>%
    as.dist() -> dist

  ps %>% sample_data() %>% data.frame() -> metadata

  ### ------

  dbRDA <- vegan::capscale(formula(paste0("dist","~",forumla)),
                           metadata,
                           add = TRUE)

  ### ------

  # overll significance of the model
  anova(dbRDA) %>%
    data.frame() -> anova_all

  # significance of different covariables
  anova(dbRDA, by = "terms") %>%
    data.frame() -> anova_terms

  # source('https://raw.githubusercontent.com/fawda123/ggord/master/R/ggord.R')

  # ggord(dbRDA, grp_in = metadata[,variables]) -> p
  autoplot(dbRDA) -> p

  ### ------

  if(!is.null(group_plot)){

    ggord(dbRDA, metadata[,group_plot],
          vec_ext = vec_ext,
          alpha = 0.5,
          ellipse_pro = 0.8,
          hull = FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) -> p2

    out <- list("plot" = p,
                "plot2"= p2,
                "dbRDA" = dbRDA,
                "anova_all" = anova_all,
                "anova_terms" = anova_terms)

  }else{

    out <- list("plot" = p,
                "dbRDA" = dbRDA,
                "anova_all" = anova_all,
                "anova_terms" = anova_terms)

  }
  ### ------

  return(out)
  ### ------

  detach("package:ggvegan", unload=TRUE);detach("package:ggord", unload=TRUE) #detach("package:plyr", unload=TRUE)

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
#' data("soilrep")
#' soilrep
#' soilrep %>%
#' phyloseq_compute_bdiv() -> test_dist
#' soilrep %>%
#' phyloseq_pairwise_dbRDA(dm = test_dist$bray, group = "Treatment")

phyloseq_pairwise_dbRDA <- function(ps,
                                    dm,
                                    RHS_formula = "Treatment"){

  require(ggvegan); require(ggord); require(BiodiversityR)

  my_dist = NULL; metadata = NULL

  as.matrix(dm)[sample_names(ps),sample_names(ps)] %>%
    as.dist() -> my_dist


  # print("this is my dist")
  # print(my_dist)
  # print("this was my dist")

  ps %>% sample_data() %>% data.frame() -> metadata

  #
  #   eval(parse(text = paste('multi_dbRDA_test <- multiconstrained(method= "capscale", my_dist ~ ', do.call(paste,
  #                                                                   c(as.list(RHS_formula_test), sep = " + ")), ",data = metadata,add = TRUE)", sep = " ")))

  multiconstrained(method="capscale", formula = formula(paste0(" my_dist ~ ", RHS_formula)),
                   data = metadata,
                   add = TRUE) -> multi_dbRDA

  out <- multi_dbRDA %>% data.frame() %>%  rownames_to_column("comp")
  return(out)

  # detach("package:ggvegan", unload=TRUE);detach("package:ggord", unload=TRUE); detach("package:BiodiversityR", unload=TRUE)
  unloadNamespace("ggvegan"); unloadNamespace("ggord"); unloadNamespace("BiodiversityR")
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
#' @examples require(phyloseq);require(tidyverse); require(vegan);sample_data(enterotype)
#' @examples enterotype %>% phyloseq::distance(method = "bray") -> bc
#' @examples sample_data(enterotype)$Var <- sample(1:40, nsamples(enterotype), replace=T)
#' @examples sample_data(enterotype)$Size <- rnorm(nsamples(enterotype), mean=176, sd=10)
#' @examples enterotype %>% phyloseq_plot_dbrda(dm = bc, grouping_column = "Project", env.variables = c("Age", "Gender", "Size", "Var"), sep = "*") -> dbrda
#' @examples

phyloseq_plot_dbrda <- function(physeq, dm, grouping_column, pvalueCutoff = 0.5, norm_method = "center_scale",
                                env.variables = NULL, num.env.variables = NULL, exclude.variables = NULL,
                                draw_species = F, nperm = 999, sep = "+")
{
  abund_table <- otu_table(physeq)
  meta_table <- data.frame(sample_data(physeq))[,c(env.variables,grouping_column)]

  complete.cases(meta_table) -> cc


  as.matrix(dm)[sample_names(physeq),sample_names(physeq)] %>%
    as.dist() -> dm

  as.matrix(dm)[cc,cc] -> dm
  meta_table[cc,] -> meta_table
  if (norm_method == "center_scale")
  {
    meta_table %>%
      mutate_if(is.numeric, scale) -> meta_table
  }
  abund_table.adonis <- vegan::adonis(formula = as.formula(paste("dm"," ~ ", paste(env.variables, collapse  = sep))),
                                      permutations = nperm,
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
  eval(parse(text = paste("sol <- vegan::capscale(dm ~ ", do.call(paste,
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

  return(list("plot"= p,
              "capscale" = sol,
              "adonis" = abund_table.adonis))
}



#' @title Plot time evolution of beta diversity metrics with respect to specified reference timepoints
#' @author Sneha Sundar and Florentin Constancias
#' @param distances a character vector of the names of the distance metrics to be plotted; valid distances should be objects in \code{bdiv_list}
#' @param bdiv_list the output of function \code{phyloseq_compute_bdiv} containing a list of computed distance matrices
#' @param physeq rarefied phyloseq object
#' @param timepoint specifies which distances should be plotted over time; can be "previous", "fixed" or "between.ref.group". More details on what these mean in notes.
#' @param group_var name of the metadata column containing group information
#' @param time_var name of the metadata column containing time information
#' @param group_to_compare name of the reference group for comparison of distances. Should be specified only if \code{timepoint} is "between.ref.group"
#' @param fixed_time the fixed time to compare distances within each group. Should be specified only if \code{timepoint} is "fixed"
#' @return a list of ggplot objects named according to the distance metric that was plotted; the plots show us how a specified kind of beta diversity metric evolves over time for each group.
#' @note
#' If \code{timepoint} is "previous",within each group in \code{group_var} ,
#' we pick out all the distances comparing current timepoint to previous
#' timepoint and plot them over time.
#'
#' If \code{timepoint} is "between.ref.group", we pick out all the distances
#' between group_to_compare and group in code{group_var} except \code{group_to_compare}
#' for which the times are the same. A simple example: if "A", "B","C" are the groups
#' and "A" is set as \code{group_to_compare}, then distances between Day.1 of B and Day.1 of A,
#' Day.2 of B and Day.2 of A,.. etc  are plotted for each group (we get one panel for "B" and another for "C").
#'
#' If \code{timepoint} is "fixed", then within each group , distances compared
#' to the \code{fixed_time} are plotted over time.
#'
#'
#'
#' @examples
#' phyloseq_plot_beta_div_wrt_timepoint(distances = c("bray","bjaccard","wjaccard"),
#'                                    bdiv_list,
#'                                    physeq=ps_polyFermS_rare,
#'                                    timepoint="between.ref.group",
#'                                    group_vartime_var="Day_from_Inoculum",
#'                                    group_to_compare="CR_UNTREATED")
#'
#'


#' @title Plot time evolution of beta diversity metrics with respect to specified reference timepoints
#' @author Sneha Sundar and Florentin Constancias
#' @param distances a character vector of the names of the distance metrics to be plotted; valid distances should be objects in \code{bdiv_list}
#' @param bdiv_list the output of function \code{phyloseq_compute_bdiv} containing a list of computed distance matrices
#' @param physeq rarefied phyloseq object
#' @param timepoint specifies which distances should be plotted over time; can be "previous", "fixed" or "between.ref.group". More details on what these mean in notes.
#' @param group_var name of the metadata column containing group information
#' @param time_var name of the metadata column containing time information
#' @param group_to_compare name of the reference group for comparison of distances. Should be specified only if \code{timepoint} is "between.ref.group"
#' @param fixed_time the fixed time to compare distances within each group. Should be specified only if \code{timepoint} is "fixed"
#' @param replicate logical indicating if there are biological replicates in the data
#' @param replicate_id_var name of the metadata column containing the replicate ids of the samples (required if \code{replicate} is TRUE)
#' @return a list of ggplot objects named according to the distance metric that was plotted; the plots show us how a specified kind of beta diversity metric evolves over time for each group.
#' @note
#' If \code{timepoint} is "previous",within each group in \code{group_var} ,
#' we pick out all the distances comparing current timepoint to previous
#' timepoint and plot them over time.
#'
#' If \code{timepoint} is "between.ref.group", we pick out all the distances
#' between group_to_compare and group in code{group_var} except \code{group_to_compare}
#' for which the times are the same. A simple example: if "A", "B","C" are the groups
#' and "A" is set as \code{group_to_compare}, then distances between Day.1 of B and Day.1 of A,
#' Day.2 of B and Day.2 of A,.. etc  are plotted for each group (we get one panel for "B" and another for "C").
#'
#' If \code{timepoint} is "fixed", then within each group , distances compared
#' to the \code{fixed_time} are plotted over time.
#'
#'
#' If \code{replicate} is TRUE, the data has biological replicates and we need to handle them
#' differently. Specify the metadata column containing replicate id information (for example if samples from 3 mice were collected, provide the metadata column specifying the mice labels )
#' We will also be plotting a boxplot rather than a time series to better visualize the distribution of distances for each day and each group.
#'
#'
#' @examples
#' phyloseq_plot_beta_div_wrt_timepoint(distances = c("bray","bjaccard","wjaccard"),
#'                                    bdiv_list,
#'                                    physeq=ps_polyFermS_rare,
#'                                    timepoint="between.ref.group",
#'                                    group_vartime_var="Day_from_Inoculum",
#'                                    group_to_compare="CR_UNTREATED",
#'                                    fixed_time=NULL,
#'                                    replicate = FALSE,
#'                                    replicate_id_var = NULL
#'                                    )


#'
#' phyloseq_plot_beta_div_wrt_timepoint(distances = c("bray","bjaccard","wjaccard"),
#'                                    bdiv_list,
#'                                    physeq=ps,
#'                                    timepoint="between.ref.group",
#'                                    group_vartime_var="Day",
#'                                    group_to_compare="H2O",
#'                                    fixed_time=NULL,
#'                                    replicate = TRUE,
#'                                    replicate_id_var = "mouse_label"
#'                                    )
#'


phyloseq_plot_beta_div_wrt_timepoint <- function(distances,
                                                 bdiv_list,
                                                 physeq,
                                                 timepoint,
                                                 group_var,
                                                 time_var,
                                                 group_to_compare=NULL,
                                                 fixed_time=NULL,
                                                 replicate = FALSE,
                                                 replicate_id_var = NULL
){


  require(phyloseq)
  require(microbiome)
  require(tidyverse)
  require(usedist)




  plot_beta_div_wrt_timepoint <- function(dist,
                                          bdiv_list,
                                          physeq,
                                          timepoint,
                                          group_var,
                                          time_var,
                                          group_to_compare=NULL,
                                          fixed_time=NULL,
                                          replicate = FALSE,
                                          replicate_id_var = NULL){

    # PREPROCESSING

    #extract distance matrix of class dist
    d.mat <- bdiv_list[[dist]]

    as.matrix(d.mat)[sample_names(physeq),sample_names(physeq)] %>%
      as.dist() -> d.mat

    #sample data as dataframe
    as(sample_data(physeq),"matrix") %>%
      data.frame(check.names=FALSE) %>%
      rownames_to_column("Sample_ID") -> sample.data

    #make time column is numeric

    if(is.numeric(sample.data[,time_var])){
      stop("Error: You need to make sure your time variable is numeric. Maybe the function parse_number() can help.")
    }


    sample.data[,time_var] <- as.numeric(sample.data[,time_var])

    #check if the labels match, if they don't throw an error
    stopifnot(all.equal(labels(d.mat), sample.data$Sample_ID))

    #extract metadata column containing group information
    item_groups <- sample.data[,group_var]

    #use library usedist to make a neat dataframe that shows for each distance in
    #`d.mat` which samples were used for calculation and the group info of both samples
    dist_df <- usedist::dist_groups(d.mat,item_groups)

    if(replicate==FALSE){
      #add the time info of the samples to `dist_df` . this is what we will use for
      #picking out the distances we want to plot
      meta_df <- sample.data %>% select(Sample_ID,.data[[time_var]])

      left_join(dist_df,
                meta_df %>%
                  dplyr::rename("varGroup1" = .data[[time_var]]
                  ),
                by = c("Item1" = "Sample_ID")) %>%
        left_join(meta_df %>%
                    dplyr::rename("varGroup2" = .data[[time_var]]),
                  by = c("Item2" = "Sample_ID")) -> dist_df


      if(timepoint=="previous")
      {
        #only within group distances are needed now
        dist_df %>%
          dplyr::filter(grepl("Within", Label)) -> dist_df


        #Create a dataframe specifying the days that we need to filter from the distance dataframe

        sample.data %>%
          select(Sample_ID,.data[[group_var]],Day1=.data[[time_var]]) %>%
          group_by(.data[[group_var]]) %>%
          arrange(Day1,.by_group=TRUE) -> days.df

        #the last day of each group should not be compared the first day of next group
        na_fills<-cumsum(days.df %>% group_size())
        #specify the reverse order so all the right comparisons are picked out
        day2<-c(days.df$Day1[-1],NA)
        day2[na_fills] <- NA

        #ungrouping adding Day2 column
        days.df <- days.df %>%
          ungroup() %>%
          mutate(Day2=day2)

        #specifying the opposite too just in case `dist_df` has it in this order
        #remmember distance between A and B is same as distance between B and A
        days.df.reverse <- days.df %>%
          rename(Day2=Day1,Day1=Day2)

        #combine both dataframes to get the complete one containing all possible day combinations of interest to us
        days.df.complete <- rbind(days.df,days.df.reverse)

        #making group info same as that of Label column in `dist_df`
        days.df.complete[,group_var] <- paste("Within",pull(days.df.complete,.data[[group_var]]))

        #Getting the right dataframe for plotting: pick out all the distances specified by `days.df.complete`
        # from dist_df
        df_plot <- semi_join(dist_df,days.df.complete,
                             by=c("Label"= group_var,"varGroup1"="Day1","varGroup2"="Day2")) %>%
          arrange(varGroup1) %>%
          arrange(varGroup2)

        #arrange and group by label
        df_plot %>%
          group_by(Label) %>%
          arrange(Label) %>%
          arrange(varGroup1,.by_group=TRUE) %>%
          arrange(varGroup2,.by_group=TRUE) -> df_plot

        #if day1 > day2 swap them to keep a consistent order: Day2 > Day1
        swap_indices<-which(df_plot$varGroup1 > df_plot$varGroup2)

        df_plot[swap_indices,c("Item1","Item2","varGroup1","varGroup2")] <- df_plot[swap_indices,c("Item2","Item1","varGroup2","varGroup1")]

        #plotting a connected scatterplot and faceting by Label
        df_plot %>%
          ggplot(aes(x = as.factor(varGroup2) , y = Distance)) +
          geom_point(size=1, alpha=0.6, aes(colour = Label, group=Label),
                     position=position_jitterdodge(dodge.width=0.9)) +
          geom_path(arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")),
                    size = 0.4, linetype = "dashed", inherit.aes = TRUE, aes(group=Label, color = Label),
                    position=position_jitterdodge(dodge.width=0.9)) +
          # geom_boxplot(aes(group = varGroup2 %>% as.factor()),
          # fill = "transparent",
          # outlier.colour = NA,alpha=0.4) +
          facet_grid(Label ~ ., scales = "fixed") +
          # ggrepel::geom_text_repel(cex=2,
          #                      aes(label= Group1),
          #                      segment.color = 'black',
          #                      segment.size = 0.5,
          #                      # nudge_x =  -4,
          #                      # nudge_y = 0,
          #                      df_all_t0 %>% filter(Day2 == 6 & metric=="wunifrac") %>% unique()) +
          theme_bw() + xlab("Day") + ylab("Distance to previous timepoint") -> plot

      }


      if(timepoint=="between.ref.group"){

        #throw error if group to compare is not specified
        if(is.null(group_to_compare)){
          stop("Error: You need to specify the name of the common group to compare. This group needs to be one of the categories in your group_var argument.")
        }

        #keeping only between sample distances with the common group. The days to be compared should be the same

        dist_df %>% dplyr::filter(!grepl("Within", Label)) ->dist_df

        dist_df %>%
          filter(Group1==group_to_compare | Group2==group_to_compare) %>%
          filter(varGroup1==varGroup2) -> dist_df


        #to keep a consistent format. reference group is group1 and the other is group2

        if(sum(dist_df$Group2==group_to_compare)!=0){
          swap_indices<-which(dist_df$Group2==group_to_compare)
          dist_df[swap_indices,c("Item1","Item2","Group1","Group2")] <- dist_df[swap_indices,c("Item2","Item1","Group2","Group1")]
          dist_df$Label <- paste("Between",dist_df$Group1,"and",dist_df$Group2)


        }

        #the dataframe we will use for plotting
        df_plot<-dist_df %>% group_by(Group2) %>% arrange(varGroup2,.by_group=TRUE)

        #plotting a connected scatterplot and faceting by Label=
        df_plot %>%
          ggplot(aes(x = as.factor(varGroup2) , y = Distance)) +
          geom_point(size=1, alpha=0.6, aes(colour = Label, group=Label),
                     position=position_jitterdodge(dodge.width=0.9)) +
          geom_path(arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")),
                    size = 0.4, linetype = "dashed", inherit.aes = TRUE, aes(group=Label, color = Label),
                    position=position_jitterdodge(dodge.width=0.9)) +
          # geom_boxplot(aes(group = varGroup2 %>% as.factor()),
          # fill = "transparent",
          # outlier.colour = NA,alpha=0.4) +
          facet_grid(Label ~ ., scales = "fixed") +
          # ggrepel::geom_text_repel(cex=2,
          #                      aes(label= Group1),
          #                      segment.color = 'black',
          #                      segment.size = 0.5,
          #                      # nudge_x =  -4,
          #                      # nudge_y = 0,
          #                      df_all_t0 %>% filter(Day2 == 6 & metric=="wunifrac") %>% unique()) +
          theme_bw() + xlab("Day") + ylab(paste("Distance to ",group_to_compare)) -> plot




      }



      if(timepoint == "fixed"){

        #throw error if fixed_time is not specified
        if(is.null(fixed_time)){
          stop("Error: You need to specify the fixed time within each group using which we will pick out distances. This time needs to be one of the times in the `time_var` column that is common to all groups.")
        }



        #only within group distances are needed now
        dist_df %>%
          dplyr::filter(grepl("Within", Label)) -> dist_df

        #get all the distances computed with a sample of `fixed_time`
        dist_df %>%
          filter(varGroup1==fixed_time | varGroup2==fixed_time) -> dist_df


        #to keep a consistent format: varGroup1 is the fixed time
        if(sum(dist_df$varGroup2==fixed_time)!=0){
          swap_indices<-which(dist_df$varGroup2==fixed_time)
          dist_df[swap_indices,c("Item1","Item2","Group1","Group2","varGroup1","varGroup2")] <-   dist_df[swap_indices,c("Item2","Item1","Group2","Group1","varGroup2","varGroup1")]

        }

        #group by label and arrange by day2
        dist_df %>% group_by(Label) %>% arrange(varGroup2,.by_group=TRUE) -> df_plot


        #plot
        df_plot %>%
          ggplot(aes(x = as.factor(varGroup2) , y = Distance)) +
          geom_point(size=1, alpha=0.6, aes(colour = Label, group=Label),
                     position=position_jitterdodge(dodge.width=0.9)) +
          geom_path(arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")),
                    size = 0.4, linetype = "dashed", inherit.aes = TRUE, aes(group=Label, color = Label),
                    position=position_jitterdodge(dodge.width=0.9)) +
          # geom_boxplot(aes(group = varGroup2 %>% as.factor()),
          # fill = "transparent",
          # outlier.colour = NA,alpha=0.4) +
          facet_grid(Label ~ ., scales = "fixed") +
          # ggrepel::geom_text_repel(cex=2,
          #                      aes(label= Group1),
          #                      segment.color = 'black',
          #                      segment.size = 0.5,
          #                      # nudge_x =  -4,
          #                      # nudge_y = 0,
          #                      df_all_t0 %>% filter(Day2 == 6 & metric=="wunifrac") %>% unique()) +
          theme_bw() + xlab("Day") + ylab(paste("Distance to",time_var,fixed_time)) -> plot


      }

      return(plot)

    }

    if(replicate == TRUE){

      #check if replicate id variable is provided . If not throw an error.
      if(is.null(replicate_id_var)){
        stop("Error: You need to specify the metadata column containing the replicate id information.")
      }


      #add the time info and replicate id info of the samples to `dist_df` . this is what we will use for
      #picking out the distances we want to plot

      meta_df <- sample.data %>% select(Sample_ID,.data[[time_var]],.data[[replicate_id_var]])

      left_join(dist_df,
                meta_df %>%
                  dplyr::rename("varGroup1" = .data[[time_var]],"rep_id1"=.data[[replicate_id_var]]
                  ),
                by = c("Item1" = "Sample_ID")) %>%
        left_join(meta_df %>%
                    dplyr::rename("varGroup2" = .data[[time_var]],"rep_id2"=.data[[replicate_id_var]]),
                  by = c("Item2" = "Sample_ID")) -> dist_df

      if(timepoint == "previous"){

        #only within group distances are needed now and we need to plot distances between the same mice
        dist_df %>%
          dplyr::filter(grepl("Within", Label)) -> dist_df
        dist_df %>% filter(rep_id1==rep_id2) -> dist_df

        sample.data %>%
          select(Sample_ID,.data[[group_var]],Day1=.data[[time_var]]) %>%
          group_by(.data[[group_var]]) %>%
          arrange(Day1,.by_group=TRUE) -> days.df

        #the last day of each group should not be compared the first day of next group
        na_fills<-cumsum(days.df %>% group_size())

        #specify the reverse order so all the right comparisons are picked out
        day2<-c(days.df$Day1[-1],NA)
        day2[na_fills] <- NA

        #ungrouping adding Day2 column
        days.df <- days.df %>%
          ungroup() %>%
          mutate(Day2=day2)

        days.df.reverse <- days.df %>%
          rename(Day2=Day1,Day1=Day2)

        #combine both dataframes to get the complete one containing all possible day combinations of interest to us
        days.df.complete <- rbind(days.df,days.df.reverse)

        #making group info same as that of Label column in `dist_df`
        days.df.complete[,group_var] <- paste("Within",pull(days.df.complete,.data[[group_var]]))

        #since we have replicate data and hence the same day comparison we need to get rid of these
        days.df.complete<-days.df.complete %>% filter(Day1!=Day2)


        #Getting the right dataframe for plotting: pick out all the distances specified by `days.df.complete`
        # from dist_df
        df_plot <- semi_join(dist_df,days.df.complete,
                             by=c("Label"= group_var,"varGroup1"="Day1","varGroup2"="Day2")) %>%
          arrange(varGroup1) %>%
          arrange(varGroup2)

        #arrange and group by label
        df_plot %>%
          group_by(Label) %>%
          arrange(Label) %>%
          arrange(varGroup1,.by_group=TRUE) %>%
          arrange(varGroup2,.by_group=TRUE) -> df_plot

        #if day1 > day2 swap them to keep a consistent order: Day2 > Day1
        swap_indices<-which(df_plot$varGroup1 > df_plot$varGroup2)

        df_plot[swap_indices,c("Item1","Item2","varGroup1","varGroup2")] <- df_plot[swap_indices,c("Item2","Item1","varGroup2","varGroup1")]



        ggplot(data=df_plot,mapping=aes(x=as.factor(varGroup2),y=Distance,color=Label)) +
          geom_boxplot(data=df_plot,mapping=aes(x=as.factor(varGroup2),y=Distance,color=Label),                 fill = NA,
                       outlier.shape = NA,
                       outlier.colour = NA) +
          geom_point(position=position_jitterdodge(jitter.width = 0.1,seed=123),aes(group=Label),size=0.1, alpha=0.4)+
          theme_bw() + xlab("Day") + ylab("Distance to previous timepoint") -> plot

      }


      if(timepoint=="between.ref.group"){

        #throw error if group to compare is not specified
        if(is.null(group_to_compare)){
          stop("Error: You need to specify the name of the common group to compare. This group needs to be one of the categories in your group_var argument.")
        }

        #keeping only between sample distances with the common group. The days to be compared should be the same

        dist_df %>% dplyr::filter(!grepl("Within", Label)) ->dist_df

        dist_df %>%
          filter(Group1==group_to_compare | Group2==group_to_compare) %>%
          filter(varGroup1==varGroup2) -> dist_df


        #to keep a consistent format. reference group is group1 and the other is group2

        if(sum(dist_df$Group2==group_to_compare)!=0){
          swap_indices<-which(dist_df$Group2==group_to_compare)
          dist_df[swap_indices,c("Item1","Item2","Group1","Group2","rep_id1","rep_id2")] <- dist_df[swap_indices,c("Item2","Item1","Group2","Group1","rep_id2","rep_id1")]
          dist_df$Label <- paste("Between",dist_df$Group1,"and",dist_df$Group2)

        }

        #the dataframe we will use for plotting
        df_plot<-dist_df %>% group_by(Group2)  %>% arrange(varGroup2,.by_group=TRUE)


        ggplot(data=df_plot,mapping=aes(x=as.factor(varGroup2),y=Distance,color=Label)) +
          geom_boxplot(data=df_plot,mapping=aes(x=as.factor(varGroup2),y=Distance,color=Label),                 fill = NA,
                       outlier.shape = NA,
                       outlier.colour = NA) +
          geom_point(position=position_jitterdodge(jitter.width = 0.1,seed=123),aes(group=Label),size=0.1, alpha=0.4) +
          theme_bw() + xlab("Day") + ylab(paste("Distance to ",group_to_compare)) -> plot
      }


      if(timepoint == "fixed"){

        #throw error if fixed_time is not specified
        if(is.null(fixed_time)){
          stop("Error: You need to specify the fixed time within each group using which we will pick out distances. This time needs to be one of the times in the `time_var` column that is common to all groups.")
        }



        #only within group distances are needed now and between the same replicates
        dist_df %>%
          dplyr::filter(grepl("Within", Label)) -> dist_df
        dist_df %>% filter(rep_id1==rep_id2) -> dist_df

        #get all the distances computed with a sample of `fixed_time`
        dist_df %>%
          filter(varGroup1==fixed_time | varGroup2==fixed_time) -> dist_df


        #to keep a consistent format: varGroup1 is the fixed time
        if(sum(dist_df$varGroup2==fixed_time)!=0){
          swap_indices<-which(dist_df$varGroup2==fixed_time)
          dist_df[swap_indices,c("Item1","Item2","Group1","Group2","varGroup1","varGroup2")] <-   dist_df[swap_indices,c("Item2","Item1","Group2","Group1","varGroup2","varGroup1")]

        }

        #group by label and arrange by day2
        dist_df %>% group_by(Label) %>% arrange(varGroup2,.by_group=TRUE) -> df_plot


        #plot
        ggplot(data=df_plot,mapping=aes(x=as.factor(varGroup2),y=Distance,color=Label)) +
          geom_boxplot(data=df_plot,mapping=aes(x=as.factor(varGroup2),y=Distance,color=Label),                 fill = NA,
                       outlier.shape = NA,
                       outlier.colour = NA) +
          geom_point(position=position_jitterdodge(jitter.width = 0.1,seed=123),aes(group=Label),size=0.1, alpha=0.4) +
          theme_bw() + xlab("Day") + ylab(paste("Distance to",time_var,fixed_time)) -> plot


      }


      return(plot)

    }

  }
  #using lapply to generate plots for a vector of distances
  res<- lapply(X=distances,FUN=plot_beta_div_wrt_timepoint,bdiv_list,
               physeq,
               timepoint,
               group_var,
               time_var,group_to_compare,
               fixed_time,
               replicate,
               replicate_id_var)

  #names of the plots will be the name of the distance metric used to compute the distance matrix
  names(res) <- distances

  return(res)
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

physeq_multi_domain_pln <- function(ps,
                                    A_filter = c("Bacteria", "Archaea"),
                                    B_filter = c("Eukaryota", "Fungi"),
                                    tax_level_filter = "Kingdom",
                                    prev_filter = 0.70,
                                    rename_taxa = FALSE,
                                    tax_glom = FALSE,
                                    formula_model = paste0("Abundance ~ 1 + offset(log(Offset))")){

  ## ------------------------------------------------------------------------
  require(tidyverse); require(phyloseq)
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using phyloseq version ", packageVersion('phyloseq'),'\n\n'))

  ## ------------------------------------------------------------------------

  if(tax_glom!=FALSE)
  {
    ps %>%
      tax_glom(taxrank = tax_glom) -> ps
  }


  if(rename_taxa){
    taxa_names(ps) <-  tax_table(ps)[,tax_glom]
  }

  ps %>%
    microbiome::core(detection=0,
                     prevalence=prev_filter) -> ps

  ## ------------------------------------------------------------------------

  #https://github.com/joey711/phyloseq/issues/1137

  # ps %>%
  #   subset_taxa(Kingdom %in% paste0(Kingdom_A_filter)) -> ps_A

  taxmat <- as(tax_table(ps), "matrix")
  taxa_A <- rownames(taxmat)[taxmat[,tax_level_filter] %in% A_filter]

  prune_taxa(taxa_A, ps) -> ps_A

  # ps %>%
  #   subset_taxa(Kingdom %in% paste0(Kingdom_B_filter)) -> ps_B

  taxa_B <- rownames(taxmat)[taxmat[,tax_level_filter] %in% B_filter]

  prune_taxa(taxa_B, ps) -> ps_B
  ## ------------------------------------------------------------------------

  ## extract counts
  counts_A <- as(phyloseq::otu_table(ps_A), "matrix") %>%
    data.frame()
  ## extract covariates (or prepare your own)
  covariates_A  <- phyloseq::sample_data(ps_A)
  ## prepare data
  my_data_A  <- prepare_data(counts = counts_A, covariates = covariates_A)

  ## extract counts
  counts_B <- as(phyloseq::otu_table(ps_B), "matrix") %>%
    data.frame()
  ## extract covariates (or prepare your own)
  covariates_B <- phyloseq::sample_data(ps_B)
  ## prepare data
  my_data_B  <- prepare_data(counts = counts_B, covariates = covariates_B)

  ## ------------------------------------------------------------------------

  don_net <- t(rbind(counts_A,counts_B))

  ## ------------------------------------------------------------------------

  offset_A <- compute_offset(t(counts_A), offset = "TSS")
  offset_B <- compute_offset(t(counts_B), offset = "TSS")

  # Offset=matrix(c(rep(offset_A, nrow(counts_A)),rep(offset_B, nrow(counts_B))),
  #               ncol = ncol(counts_A) + ncol(counts_B), nrow = nrow(counts_A) + nrow(counts_B), byrow=TRUE)
  #
  Offset = matrix(c(rep(offset_A, nrow(counts_A)),rep(offset_B, nrow(counts_B))),
                  dim(don_net)[1], nrow=dim(don_net)[2], byrow=TRUE)

  Offset=t(Offset)
  ## ------------------------------------------------------------------------

  data_cov <- prepare_data(counts = don_net, covariates = covariates_B)

  data_cov$Offset <- Offset
  ## ------------------------------------------------------------------------

  models_net <- PLNnetwork(as.formula(formula_model), data = data_cov)

  ## ------------------------------------------------------------------------

  models_net %>%
    plot("diagnostic") -> diag_p

  ## ------------------------------------------------------------------------

  models_net %>%
    plot() -> model_p

  ## ------------------------------------------------------------------------

  models_net %>%
    coefficient_path(corr = TRUE) %>%
    ggplot(aes(x = Penalty, y = Coeff, group = Edge, colour = Edge)) +
    geom_line(show.legend = FALSE) +  coord_trans(x="log10") + theme_bw() -> model_p2


  ## ------------------------------------------------------------------------


  models_net %>%
    getBestModel("StARS") -> model_StARS # if StARS is requested, stabiltiy selection is performed if needed

  models_net %>%
    getBestModel("BIC") -> model_BIC# if StARS is requested, stabiltiy selection is performed if needed

  ## ------------------------------------------------------------------------


  ## ------------------------------------------------------------------------

  out <- list("models" = models_net,
              "diag_p" = diag_p,
              "model_p" = model_p,
              "model_p2" = model_p2,
              "model_StARS" = model_StARS,
              "model_BIC" = model_BIC,
              "lambda" = models_net$penalties,
              "my_graph" = plot(model_StARS, plot = FALSE),
              "model_p3" = plot(model_StARS),
              "model_p4" = data.frame(
                fitted   = as.vector(fitted(model_StARS)),
                observed = as.vector(data_cov$Abundance)
              ) %>%
                ggplot(aes(x = observed, y = fitted)) +
                geom_point(size = .5, alpha =.25 ) +
                scale_x_log10(limits = c(1,1000)) +
                scale_y_log10(limits = c(1,1000)) +
                theme_bw() + annotation_logticks())

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
#'

in_vitro_mIMT_STABvsTreat <- function(physeq,
                                      var = "treatment",
                                      dlist,
                                      group,
                                      formula,
                                      strata,
                                      group_color,
                                      group_shape,
                                      group_alpha,
                                      m = "PCoA")
{

  ps <- prune_samples(get_variable(physeq, var) == group,
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

  p = p + facet_wrap(distance ~ ., scales="free")

  p = p + ggtitle(paste0(m," using various distance metrics ", "- ", group, " Group")) +
    theme_light() #+ scale_color_viridis_d()

  p + scale_color_viridis_d() +
    scale_alpha_continuous(range = c(0.6, 1),
                           breaks = c(20,30,40)) +
    scale_shape_manual(values = c(1,19,0,15))

  output = list("plot" = p,
                "PERMANOVA" = adonis_tmp)

  return(output)
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

phyloseq_generate_pcoa_per_variables <- function(tmp,
                                                 group,
                                                 dist){


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
      filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> filter_tmp

    filter_tmp %>%
      phyloseq_plot_bdiv(dlist = dist,
                         seed = 123,
                         axis1 = 1,
                         axis2 = 2) -> out[[tp]]
  }

  return(out)
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note todo:     filter(rlang::eval_tidy(rlang::parse_expr(filtering_expr))) -> sample_pw_meta
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#' A column 'Stab_Treat' in your metadata should describe if samples were collected during the Stab or Treatment phase
#' A column 'Period' in your metadata should describe the treatment Period / fermentation / repetition using same donor.
#' A column 'Treatment' in your metadata should describe the Treatment applied.
#' Only pairwise ratio will be computed for each ASV (or taxa if you specify tax_rank) within the same Treatment factors (Treatment metadata), within the same Period (Period metadata) and between TREAT vs STAB (Stab_Treat metadata, STAB is the reference)


plot_ratio_Stab_Treat_taxa <- function(ps_new,
                                       transform = "compositional", # percentage 0,1
                                       # tax_rank = "Genus",
                                       taxa_sel = NULL, #to be implemented
                                       detection = 0, # Detection threshold for absence/presence
                                       prevalence = 0.20, # Prevalence threshold (in [0, 1]).
                                       n_filter = 4, # Minimum number of ratio per ASV/Taxa to be kept : to plot - perform stats
                                       meta_sel = c("sample_name","Stab_Treat", "Period", "Reactor","Timepoint", "Treatment","Period"), # selection of metadata to be exported with the data
                                       plot = TRUE,
                                       stats = FALSE,
                                       ref_group_stat = "NA",
                                       to_remove = c("unknown", "Incertae Sedis")){

  #### ------------- get pairwise sample combinations

  ps_new %>%
    sample_names() %>%
    as.vector() -> samples_names

  tibble::tibble(Sample_A = samples_names,
                 Sample_B = samples_names) %>%
    tidyr::expand(Sample_A, Sample_B) %>%
    dplyr::filter(Sample_A != Sample_B) -> sample_pw  # remove self comparaisons

  #### ------------- transformation prevalence filtering - if no taconomic agglomeration

  if(is.null(tax_rank)){

    ps_new %>%
      microbiome::transform(transform = transform) %>%
      microbiome::core(., detection = detection, prevalence = prevalence ) -> ps_new

  }
  #### ------------- Taxonomic agglomeration and transformation prevalence filtering

  if(!is.null(tax_rank)){
    ps_new %>%
      speedyseq::tax_glom(taxrank = tax_rank) %>%
      microbiome::transform(transform = transform) %>%
      microbiome::core(., detection = detection, prevalence = prevalence ) -> ps_new

    prune_taxa(data.frame(tax_table(ps_new)[,tax_rank])
               %>%  dplyr::filter(!get(tax_rank) %in% to_remove) %>% rownames(),ps_new) -> ps_new

    taxa_names(ps_new) <-  tax_table(ps_new)[,tax_rank]

  }

  #### ------------- Extracting metadata based on meta_sel parameter

  ps_new %>%
    sample_data() %>%
    data.frame() %>%
    dplyr::select(one_of(meta_sel)) -> meta_data

  meta_data_A <- meta_data
  meta_data_B <- meta_data

  names(meta_data_A) <- paste0(names(meta_data), "_A")
  names(meta_data_B) <- paste0(names(meta_data), "_B")

  #### ------------- Joining metadata with pairwise sample comparaisons

  sample_pw %>%
    left_join(meta_data_A,
              by = c("Sample_A" = "sample_name_A")) %>%
    left_join(meta_data_B,
              by = c("Sample_B" = "sample_name_B")) %>%
    filter(Stab_Treat_A  %in% "STAB" & # Sample_A has to be STAB  (Stab_Treat_A colum)
             Stab_Treat_B == "TREAT" & # and Sample_B has to be TREAT (Stab_Treat_B colum)
             Period_A == Period_B, # and Period_A (of Sample_A has to be the same as the period of Sample_B - no between period ratios)
           Treatment_A == Treatment_B) -> sample_pw_meta # only within treatment ratios

  #### ------------- Exporting the OTU 'abundance' data

  ps_new %>%
    speedyseq::psmelt() %>%
    dplyr::select(Sample, OTU, Abundance) -> data

  data_A <- data
  data_B <- data

  names(data_A) <- paste0(names(data), "_A")
  names(data_B) <- paste0(names(data), "_B")

  #### ------------- Joining OTU abundance with pairwise sample comparisons + metadata already filter above

  sample_pw_meta %>%
    left_join(.,
              data_A,
              by = c("Sample_A" = "Sample_A")) %>%
    left_join(.,
              data_B,
              by = c("Sample_B" = "Sample_B")) %>%
    dplyr::filter(OTU_A == OTU_B) %>% # only within 'OTU' ratio - could be Genus...
    dplyr::mutate(ratio = Abundance_B  / Abundance_A ) %>%  # computing ratio
    dplyr::mutate(log10_ratio = log10(ratio)) -> sample_pw_meta_tax

  #### -------------

  sample_pw_meta_tax %>%
    dplyr::filter(is.finite(log10_ratio)) %>%
    group_by(Period_A, Reactor_A, Treatment_A, OTU_A) %>%
    add_count() %>% # select(n, OTU_A)
    filter(n() >= n_filter) %>% # keeping only OTU (Genus, ...) with at least n_filter pairwise comparaisons
    ungroup() -> sample_pw_meta_tax_filt

  #### -------- plot
  if(plot == TRUE){

    boxplot_ratio <- sample_pw_meta_tax_filt %>%
      ggplot(., aes(x = OTU_A, y = log10(ratio))) +
      geom_boxplot(aes(color = Treatment_A, fill = Treatment_A),
                   outlier.shape = NA,
                   outlier.colour = NA,
                   alpha = 0.4, position = position_dodge2(preserve = "single")) +
      geom_point(aes(color = Treatment_A, shape = Period_A),
                 alpha = 0.3, size  = 1.25,  position = position_jitterdodge(dodge.width = 0.8)) +
      # facet_grid(. ~ OTU_A  , scales = "free") +
      ylab(paste0("log10 ratio Treat/Stab")) + xlab(NULL) + #ylim(c(0,1)) +
      theme_light() + #theme(legtheme(legend.position = "none") +
      ggpubr::rotate_x_text(45)

    out <- list("df" = sample_pw_meta_tax,
                "df_filtered" = sample_pw_meta_tax_filt,
                "plot" = boxplot_ratio)
  }

  #### -------- Stats
  if(stats == TRUE){

    sample_pw_meta_tax_filt %>%
      group_by(OTU_A, Period_A) %>%
      rstatix::wilcox_test(log10_ratio ~ Treatment_A,
                           data = . ) %>%  #,
      # ref.group = "Control") %>%
      filter(group1 == !!ref_group_stat) %>%
      rstatix::adjust_pvalue(method = "fdr") %>%
      rstatix::add_significance("p.adj") -> boxplot_ratio_stats

    out <- list("df" = sample_pw_meta_tax,
                "df_filtered" = sample_pw_meta_tax_filt,
                "plot" = boxplot_ratio,
                "stats" = stats)
  }
  if(stats == FALSE & plot == FALSE){
    out <- list("df" = sample_pw_meta_tax,
                "df_filtered" = sample_pw_meta_tax_filt)
  }


  return(out)
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note todo:     filter(rlang::eval_tidy(rlang::parse_expr(filtering_expr))) -> sample_pw_meta
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'

plot_ratio_Stab_Treat_meta <- function(ps_new,
                                       meta_sel = c("sample_name","Stab_Treat", "Period", "Reactor","Timepoint", "Treatment","Period"),
                                       meta_plot = "Clostridium.perfringens",
                                       # plot = TRUE,
                                       # stats = FALSE,
                                       ref_group_stat = "NA"){

  #### ------------- get pairwise sample combinations

  ps_new %>%
    sample_names() %>%
    as.vector() -> samples_names

  tibble::tibble(Sample_A = samples_names,
                 Sample_B = samples_names) %>%
    tidyr::expand(Sample_A, Sample_B) %>%
    dplyr::filter(Sample_A != Sample_B) -> sample_pw  # remove self comparaisons

  #### ------------- Extracting metadata based on meta_sel parameter

  ps_new %>%
    sample_data() %>%
    data.frame() %>%
    dplyr::select(one_of(meta_sel, meta_plot)) %>%
    dplyr::rename(meta_plot = all_of(meta_plot)) -> meta_data

  meta_data_A <- meta_data
  meta_data_B <- meta_data

  names(meta_data_A) <- paste0(names(meta_data), "_A")
  names(meta_data_B) <- paste0(names(meta_data), "_B")

  #### ------------- Joining metadata with pairwise sample comparaisons

  sample_pw %>%
    left_join(meta_data_A,
              by = c("Sample_A" = "sample_name_A")) %>%
    left_join(meta_data_B,
              by = c("Sample_B" = "sample_name_B")) %>%
    filter(Stab_Treat_A  %in% "STAB" & # Sample_A has to be STAB  (Stab_Treat_A colum)
             Stab_Treat_B == "TREAT" & # and Sample_B has to be TREAT (Stab_Treat_B colum)
             Period_A == Period_B, # and Period_A (of Sample_A has to be the same as the period of Sample_B - no between period ratios)
           Treatment_A == Treatment_B) -> sample_pw_meta # only within treatment ratios

  #### ------------- Joining OTU abundance with pairwise sample comparisons + metadata already filter above
  sample_pw_meta %>%
    dplyr::mutate(ratio = meta_plot_A / meta_plot_B) %>%  # computing ratio
    dplyr::mutate(log10_ratio = log10(ratio)) -> sample_pw_meta

  #### -------------

  # sample_pw_meta %>%
  #   dplyr::filter(is.finite(log10_ratio)) %>%
  #   group_by(Period_A, Reactor_A, Treatment_A, meta_plot_A) %>%
  #   add_count() %>% # select(n, OTU_A)
  #   filter(n() >= n_filter) %>% # keeping only OTU (Genus, ...) with at least n_filter pairwise comparaisons
  #   ungroup() -> sample_pw_meta_filt

  #### -------- plot
  if(plot == TRUE){

    boxplot_ratio <- sample_pw_meta %>%
      ggplot(., aes(x = Treatment_A, y = log10(ratio))) +
      geom_boxplot(aes(color = Treatment_A, fill = Treatment_A),
                   outlier.shape = NA,
                   outlier.colour = NA,
                   alpha = 0.4, position = position_dodge2(preserve = "single")) +
      geom_point(aes(color = Treatment_A, shape = Period_A),
                 alpha = 0.3, size  = 1.25,  position = position_jitterdodge(dodge.width = 0.8)) +
      # facet_grid(. ~ OTU_A  , scales = "free") +
      ylab(paste0("log10 ratio Treat/Stab ", meta_plot)) + xlab(NULL) + #ylim(c(0,1)) +
      theme_light() + #theme(legtheme(legend.position = "none") +
      ggpubr::rotate_x_text(45)

    out <- list("df" = sample_pw_meta,
                "plot" = boxplot_ratio)
  }

  #### -------- Stats
  if(stats == TRUE){

    sample_pw_meta %>%
      group_by(meta_plot_A, Period_A) %>%
      rstatix::wilcox_test(log10_ratio ~ Treatment_A,
                           data = . ) %>%  #,
      # ref.group = "Control") %>%
      filter(group1 == !!ref_group_stat) %>%
      rstatix::adjust_pvalue(method = "fdr") %>%
      rstatix::add_significance("p.adj") -> boxplot_ratio_stats

    out <- list("df" = sample_pw_meta,
                "plot" = boxplot_ratio,
                "stats" = stats)
  }
  if(stats == FALSE & plot == FALSE){
    out <- list("df" = sample_pw_meta)
  }


  return(out)
}



#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note https://github.com/USFOneHealthCodeathon2020/Team1_MicroPowerPlus/blob/master/PERMANOVA.Rmd
#' @note https://htmlpreview.github.io/?https://github.com/USFOneHealthCodeathon2020/Team1_MicroPowerPlus/blob/master/PERMANOVA.html
#' @note Expect the terms to be under the terms column and not as rownames.
#' @return .
#' @export
#' @examples
#'library(phyloseq)
#'data(enterotype)
#'phyloseq_compute_bdiv(enterotype) -> dist
#'enterotype %>%
#'subset_samples(!is.na(Gender)) %>%
#'phyloseq_adonis2(dm = dist$bray, physeq = ., formula = "Gender") -> adonis2
#'adonis2 %>% adonis_OmegaSq(partial = FALSE)

adonis_OmegaSq <- function(aov_tab, partial = TRUE){

  # if(is.numeric(colSums(aov_tab, na.rm = TRUE)))
  # aov_tab %>%
  #   rownames_to_column('terms') -> aov_tmp

  ####---------------------- Compute MeanSqs

  aov_tab %>%
    mutate(MeanSqs = SumOfSqs / Df) -> aov_tmp

  ####---------------------- Identify MS_res SS_tot and N

  aov_tmp %>%
    filter(.[[1]]  == "Residual") %>%
    pull(MeanSqs) -> MS_res

  aov_tmp %>%
    filter(.[[1]]  == "Total") %>%
    pull(SumOfSqs) -> SS_tot

  aov_tmp %>%
    filter(.[[1]]  == "Total") %>%
    pull(Df) + 1 -> N

  ####---------------------- Run (partial) Omega Square

  if(partial == TRUE){
    omega <- apply(aov_tmp %>% column_to_rownames('terms'), 1, function(x) (x["Df"]*(x["MeanSqs"]-MS_res))/(x["Df"]*x["MeanSqs"]+(N-x["Df"])*MS_res))
    aov_tmp$parOmegaSq <- c(omega[1:(length(omega)-2)], NA, NA)
    cn_order <- c("Df", "SumOfSqs", "MeanSqs", "F", "R2", "parOmegaSq", "Pr(>F)")
  } else {
    omega <- apply(aov_tmp %>% column_to_rownames('terms'), 1, function(x) (x["SumsOfSqs"]-x["Df"]*MS_res)/(SS_tot+MS_res))
    aov_tmp$OmegaSq <- c(omega[1:(length(omega)-2)], NA, NA)
    cn_order <- c("Df", "SumOfSqs", "MeanSqs", "F", "R2", "OmegaSq", "Pr(>F)")
  }

  ####---------------------- Reorder the table and return the output

  aov_tmp %>%
    select(terms, one_of(cn_order)) -> aov_tmp

  return(aov_tmp)
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note mntd.obs.z 	 Standardized effect size of MNTD vs. null communities (= (mntd.obs - mntd.rand.mean) / mntd.rand.sd, equivalent to -NTI)
#' @return .
#' @export
#' @examples
#' library(phyloseq);library(tidyverse);data("GlobalPatterns")
#'
#' system.time(GlobalPatterns %>% subset_samples(SampleType == "Feces") %>% filter_taxa(function(x) sum(x > 0) > 2, TRUE) %>% phyloseq_comdist_parallel(., cores = 6) -> out_phyloseq_comdistNTI_parallel)

#' https://rfunctions.blogspot.com/2012/07/standardized-effect-size-nearest.html


phyloseq_comdist_parallel <- function(physeq, abundance.weighted = FALSE, exclude.conspecifics = FALSE, cores = 1, progress = TRUE){

  ####---------------------- Load R package
  require(picante); require(doSNOW)

  ####---------------------- Extract data
  physeq %>%
    microbiome::transform(transform = "compositional") -> physeq

  as(otu_table(physeq), "matrix") %>%
    t() %>%
    as.data.frame() -> comm


  physeq@phy_tree %>%
    cophenetic() -> dis

  ####----------------------

  dat <- match.comm.dist(comm, dis)
  comm <- dat$comm
  dis <- dat$dist
  N <- dim(comm)[1]

  sppInSamples <- apply(comm,1,function(x) names(which(x > 0)))

  if(progress){
    pb <- txtProgressBar(max = (N-1), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  } else {
    opts <- NULL
  }

  if(cores == 1) {
    registerDoSEQ() }
  else {
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  }

  i <- NULL
  comdisnt <- foreach (i = 1:(N-1),.combine = rbind, .options.snow = opts) %dopar% {

    comdisnt.sub <- as.numeric(rep(NA,N))

    for (j in (i + 1):N) {

      sppInSample1 <- sppInSamples[[i]]
      sppInSample2 <- sppInSamples[[j]]

      if ((length(sppInSample1) >= 1) && (length(sppInSample2) >= 1)) {
        sample.dis <- dis[sppInSample1, sppInSample2, drop = FALSE]
        if (exclude.conspecifics) {
          sample.dis[sample.dis == 0] <- NA
        }
        sample1NT <- apply(sample.dis, 1, min, na.rm = TRUE)
        sample1NT[sample1NT == Inf] <- NA
        sample2NT <- apply(sample.dis, 2, min, na.rm = TRUE)
        sample2NT[sample2NT == Inf] <- NA

        if (abundance.weighted) {
          sample1.weights <- as.numeric(comm[i, sppInSample1])
          sample2.weights <- as.numeric(comm[j, sppInSample2])
          if (any(is.na(sample1NT))) {
            miss <- which(is.na(sample1NT))
            sample1NT <- sample1NT[-miss]
            sample1.weights <- sample1.weights[-miss]
            sample1.weights <- sample1.weights/sum(sample1.weights)
          }
          if (any(is.na(sample2NT))) {
            miss <- which(is.na(sample2NT))
            sample2NT <- sample2NT[-miss]
            sample2.weights <- sample2.weights[-miss]
            sample2.weights <- sample2.weights/sum(sample2.weights)
          }
          sampleNT <- c(sample1NT, sample2NT)
          sample.weights <- c(sample1.weights, sample2.weights)
          comdisnt.sub[j] <- weighted.mean(sampleNT, sample.weights,  na.rm = TRUE)
        }
        else {
          comdisnt.sub[j] <- mean(c(sample1NT, sample2NT), na.rm = TRUE)
        }
      }
      else {
        comdisnt.sub[j] <- NA
      }
    }
    return(comdisnt.sub)
  }
  if(cores != 1) stopCluster(cl)
  comdisnt <- rbind(comdisnt,rep(NA,N))

  rownames(comdisnt) <- colnames(comdisnt) <- rownames(comm)
  return(as.dist(t(comdisnt)))

  detach("package:picante", unload=TRUE); detach("package:doSNOW", unload=TRUE)

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
#' library(phyloseq);library(tidyverse);data("GlobalPatterns")
#'
#' phyloseq_comdist_SES_NTI_parallel(GlobalPatterns %>%  subset_samples(SampleType == "Skin") %>% filter_taxa(function(x) sum(x > 0) > 2, TRUE), cores = 4) -> out
#' https://rfunctions.blogspot.com/2012/07/standardized-effect-size-nearest.html


phyloseq_comdist2_SES_NTI_parallel <- function(physeq, method = "swap", fixedmar = "both", shuffle = "both", strata = NULL, mtype = "count", burnin = 0, thin = 1,
                                               abundance.weighted = FALSE, exclude.conspecifics = FALSE, runs = 999, cores = 1){


  ####---------------------- Load R package
  require(picante); require(vegan)

  ####---------------------- Extract data
  # physeq %>%
  #   microbiome::transform(transform = "compositional") -> physeq

  as(otu_table(physeq), "matrix") %>%
    t() %>%
    as.data.frame() -> samp


  physeq@phy_tree %>%
    cophenetic() %>%
    as.matrix() -> dis

  ####----------------------

  comdistnt.obs <- as.matrix(comdistnt.par(samp, dis, abundance.weighted = abundance.weighted, exclude.conspecifics = exclude.conspecifics, cores = cores, progress = FALSE))

  if(is.null(method)) {
    comdistnt.rand <- replicate(runs, as.matrix(comdist_nti_par(permatfull(samp, fixedmar = fixedmar, shuffle = shuffle, strata = strata, mtype = mtype, times = 1)$perm[[1]], dis, abundance.weighted, exclude.conspecifics, cores = cores, progress = FALSE)), simplify = FALSE)
  } else {
    comdistnt.rand <- replicate(runs, as.matrix(comdist_nti_par(permatswap(samp, method = method, fixedmar = fixedmar, shuffle = shuffle, strata = strata, mtype = mtype, burnin = burnin, thin = thin, times = 1)$perm[[1]], dis, abundance.weighted, exclude.conspecifics, cores = cores, progress = FALSE)), simplify = FALSE)
  }

  comdistnt.rand.mean <- apply(X = simplify2array(comdistnt.rand), MARGIN = 1:2, FUN = mean, na.rm = TRUE)

  comdistnt.rand.sd <- apply(X = simplify2array(comdistnt.rand), MARGIN = 1:2, FUN = sd, na.rm = TRUE)

  comdistnt.obs.z <- (comdistnt.obs - comdistnt.rand.mean)/comdistnt.rand.sd

  comdistnt.obs.rank <- apply(X = simplify2array(c(list(comdistnt.obs),comdistnt.rand)), MARGIN = 1:2, FUN = rank)[1,,]
  comdistnt.obs.rank <- ifelse(is.na(comdistnt.rand.mean), NA, comdistnt.obs.rank)
  diag(comdistnt.obs.rank) <- NA

  comdistnt.obs.p <- comdistnt.obs.rank/(runs + 1)

  list(ntaxa = specnumber(samp), comdistnt.obs = comdistnt.obs, comdistnt.rand.mean = comdistnt.rand.mean,
       comdistnt.rand.sd = comdistnt.rand.sd, comdistnt.obs.rank = comdistnt.obs.rank, comdistnt.obs.z = comdistnt.obs.z, comdistnt.obs.p = comdistnt.obs.p, runs = runs)
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note mntd.obs.z 	 Standardized effect size of MNTD vs. null communities (= (mntd.obs - mntd.rand.mean) / mntd.rand.sd, equivalent to -NTI)
#' @note TO CHECK: https://search.r-project.org/CRAN/refmans/iCAMP/html/bNTIn.p.html & RC.pc with microeco::
#' @note ... & https://docs.ropensci.org/phylocomr/reference/ph_comdist.html
#' @return
#' @export
#' @examples


#' MicEco::ses.comdistnt()
#' library(phyloseq);library(tidyverse);data("GlobalPatterns")
#'
#'system.time(GlobalPatterns %>% subset_samples(SampleType == "Feces") %>% filter_taxa(function(x) sum(x > 0) > 2, TRUE) %>% phyloseq_bNTI_parallel(.,   abundance_weighted = TRUE, exclude_conspecifics = FALSE, null.model = "taxa.labels", cores = 4, runs = 4, iterations = 4) -> out_phyloseq_bNTI_parallel)

phyloseq_bNTI_parallel <- function(physeq, null.model = c("taxa.labels", "richness",
                                                          "frequency", "sample.pool", "phylogeny.pool", "independentswap",
                                                          "trialswap"), abundance_weighted = TRUE, exclude_conspecifics = FALSE,
                                   runs = 999, iterations = 1000, cores = 1)
{


  ####---------------------- Load R package
  require(picante); require(vegan);require(doSNOW)

  ####---------------------- Extract data
  # physeq %>%
  #   microbiome::transform(transform = "compositional") -> physeq

  as(otu_table(physeq), "matrix") %>%
    t() %>%
    as.data.frame() -> samp


  physeq@phy_tree %>%
    cophenetic() %>%
    as.matrix() -> dis

  ####----------------------

  comdistnt.obs <- as.matrix(comdist_nti_par(samp, dis, abundance_weighted = abundance_weighted,
                                             exclude_conspecifics = exclude_conspecifics, cores = cores,
                                             progress = FALSE))
  null.model <- match.arg(null.model)
  comdistnt.rand <- switch(null.model, taxa.labels = replicate(runs,
                                                               as.matrix(comdist_nti_par(samp, taxaShuffle(dis), abundance_weighted = abundance_weighted,
                                                                                         exclude_conspecifics = exclude_conspecifics, cores = cores,
                                                                                         progress = FALSE)), simplify = FALSE), richness = replicate(runs,
                                                                                                                                                     as.matrix(comdist_nti_par(randomizeMatrix(samp, null.model = "richness"),
                                                                                                                                                                               dis, abundance_weighted, exclude_conspecifics = exclude_conspecifics,
                                                                                                                                                                               cores = cores, progress = FALSE)), simplify = FALSE),
                           frequency = replicate(runs, as.matrix(comdist_nti_par(randomizeMatrix(samp,
                                                                                                 null.model = "frequency"), dis, abundance_weighted,
                                                                                 exclude_conspecifics = exclude_conspecifics, cores = cores,
                                                                                 progress = FALSE)), simplify = FALSE), sample.pool = replicate(runs,
                                                                                                                                                as.matrix(comdist_nti_par(randomizeMatrix(samp, null.model = "richness"),
                                                                                                                                                                          dis, abundance_weighted, exclude_conspecifics = exclude_conspecifics,
                                                                                                                                                                          cores = cores, progress = FALSE)), simplify = FALSE),
                           phylogeny.pool = replicate(runs, as.matrix(comdist_nti_par(randomizeMatrix(samp,
                                                                                                      null.model = "richness"), taxaShuffle(dis), abundance_weighted,
                                                                                      exclude_conspecifics = exclude_conspecifics, cores = cores,
                                                                                      progress = FALSE)), simplify = FALSE), independentswap = replicate(runs,
                                                                                                                                                         as.matrix(comdist_nti_par(randomizeMatrix(samp, null.model = "independentswap",
                                                                                                                                                                                                   iterations), dis, abundance_weighted, exclude_conspecifics = exclude_conspecifics,
                                                                                                                                                                                   cores = cores, progress = FALSE)), simplify = FALSE),
                           trialswap = replicate(runs, as.matrix(comdist_nti_par(randomizeMatrix(samp,
                                                                                                 null.model = "trialswap", iterations), dis, abundance_weighted,
                                                                                 exclude_conspecifics = exclude_conspecifics, cores = cores,
                                                                                 progress = FALSE)), simplify = FALSE))
  comdistnt.rand.mean <- apply(X = simplify2array(comdistnt.rand),
                               MARGIN = 1:2, FUN = mean, na.rm = TRUE)
  comdistnt.rand.sd <- apply(X = simplify2array(comdistnt.rand),
                             MARGIN = 1:2, FUN = sd, na.rm = TRUE)
  comdistnt.obs.z <- (comdistnt.obs - comdistnt.rand.mean)/comdistnt.rand.sd
  comdistnt.obs.rank <- apply(X = simplify2array(c(list(comdistnt.obs),
                                                   comdistnt.rand)), MARGIN = 1:2, FUN = rank)[1, , ]
  comdistnt.obs.rank <- ifelse(is.na(comdistnt.rand.mean),
                               NA, comdistnt.obs.rank)
  diag(comdistnt.obs.rank) <- NA
  comdistnt.obs.p <- comdistnt.obs.rank/(runs + 1)

  # bNTI <-   weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals)


  out <- list(ntaxa = specnumber(samp), comdistnt_obs = comdistnt.obs,
              comdistnt_rand.mean = comdistnt.rand.mean, comdistnt_rand_sd = comdistnt.rand.sd,
              comdistnt_obs.rank = comdistnt.obs.rank, comdistnt_obs_z = comdistnt.obs.z,
              comdistnt_obs.p = comdistnt.obs.p, bNTI = -1 * comdistnt.obs.z, runs = runs)


  return(out)

  detach("package:doSNOW", unload=TRUE); detach("package:picante", unload=TRUE)

}


comdist_nti_par <- function (comm, dis, abundance_weighted = FALSE, exclude_conspecifics = FALSE,
                             cores = 1, progress = TRUE)
{
  dat <- match.comm.dist(comm, dis)
  comm <- dat$comm
  dis <- dat$dist
  N <- dim(comm)[1]
  comm <- decostand(comm, method = "total", MARGIN = 1)
  sppInSamples <- apply(comm, 1, function(x) names(which(x >
                                                           0)))
  if (progress) {
    pb <- txtProgressBar(max = (N - 1), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }
  else {
    opts <- NULL
  }
  if (cores == 1) {
    registerDoSEQ()
  }
  else {
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  }
  i <- NULL
  comdisnt <- foreach(i = 1:(N - 1), .combine = rbind, .options.snow = opts) %dopar%
    {
      comdisnt.sub <- as.numeric(rep(NA, N))
      for (j in (i + 1):N) {
        sppInSample1 <- sppInSamples[[i]]
        sppInSample2 <- sppInSamples[[j]]
        if ((length(sppInSample1) >= 1) && (length(sppInSample2) >=
                                            1)) {
          sample.dis <- dis[sppInSample1, sppInSample2,
                            drop = FALSE]
          if (exclude_conspecifics) {
            sample.dis[sample.dis == 0] <- NA
          }
          sample1NT <- apply(sample.dis, 1, min, na.rm = TRUE)
          sample1NT[sample1NT == Inf] <- NA
          sample2NT <- apply(sample.dis, 2, min, na.rm = TRUE)
          sample2NT[sample2NT == Inf] <- NA
          if (abundance_weighted) {
            sample1.weights <- as.numeric(comm[i, sppInSample1])
            sample2.weights <- as.numeric(comm[j, sppInSample2])
            if (any(is.na(sample1NT))) {
              miss <- which(is.na(sample1NT))
              sample1NT <- sample1NT[-miss]
              sample1.weights <- sample1.weights[-miss]
              sample1.weights <- sample1.weights/sum(sample1.weights)
            }
            if (any(is.na(sample2NT))) {
              miss <- which(is.na(sample2NT))
              sample2NT <- sample2NT[-miss]
              sample2.weights <- sample2.weights[-miss]
              sample2.weights <- sample2.weights/sum(sample2.weights)
            }
            sampleNT <- c(sample1NT, sample2NT)
            sample.weights <- c(sample1.weights, sample2.weights)
            comdisnt.sub[j] <- weighted.mean(sampleNT,
                                             sample.weights, na.rm = TRUE)
          }
          else {
            comdisnt.sub[j] <- mean(c(sample1NT, sample2NT),
                                    na.rm = TRUE)
          }
        }
        else {
          comdisnt.sub[j] <- NA
        }
      }
      return(comdisnt.sub)
    }
  if (cores != 1)
    stopCluster(cl)
  comdisnt <- rbind(comdisnt, rep(NA, N))
  rownames(comdisnt) <- colnames(comdisnt) <- rownames(comm)
  return(as.dist(t(comdisnt)))
}

#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note mntd.obs.z 	 Standardized effect size of MNTD vs. null communities (= (mntd.obs - mntd.rand.mean) / mntd.rand.sd, equivalent to -NTI)
#' @note TO CHECK: https://search.r-project.org/CRAN/refmans/iCAMP/html/bNTIn.p.html & RC.pc with microeco::
#' @note ... & https://docs.ropensci.org/phylocomr/reference/ph_comdist.html
#' @return .
#' @export
#' @examples
#'
#'library(phyloseq);library(tidyverse);data("GlobalPatterns")
#'
#'GlobalPatterns %>% subset_samples(SampleType == "Feces") %>% filter_taxa(function(x) sum(x > 0) > 2, TRUE) %>% phyloseq_iCAMP_bNTIn_par(physeq = ., rand = 4, weighted = TRUE, exclude.consp = FALSE, nworker = 6, detail.null= FALSE) -> out_phyloseq_iCAMP_bNTIn_par
#'
#'require(iCAMP)
#'data("example.data")
#'comm=example.data$comm
#'pd=example.data$pd
#'nworker=2 # parallel computing thread number
#'rand.time=4 # usually use 1000 for real data.
#'bNTI=bNTIn.p(comm=comm, dis=pd, nworker = nworker, memo.size.GB = 50,
#'             weighted = TRUE, exclude.consp = FALSE, rand = rand.time,
#'             output.bMNTD = FALSE, sig.index = "SES", unit.sum = NULL,
#'             correct.special = TRUE, detail.null = FALSE,
#'             special.method = "MNTD")

phyloseq_iCAMP_bNTIn_par <- function(physeq, nworker = 4, memo.size.GB = 50,
                                     weighted = TRUE, exclude.consp = FALSE,
                                     rand = 1000, output.bMNTD = FALSE,
                                     sig.index=c("SES"),
                                     unit.sum = NULL, correct.special = FALSE,
                                     detail.null=TRUE,
                                     special.method=c("MNTD"),
                                     ses.cut=1.96,rc.cut=0.95,conf.cut=0.975,
                                     dirichlet = FALSE){


  ####---------------------- Load R package
  require(iCAMP); require(vegan)

  ####---------------------- Extract data
  # physeq %>%
  #   microbiome::transform(transform = "compositional") -> physeq

  as(otu_table(physeq), "matrix") %>%
    t() %>%
    as.matrix() -> commps


  physeq@phy_tree %>%
    cophenetic() %>%
    as.matrix() -> dis

  ####----------------------


  bNTIn.p(comm = commps, dis = dis, nworker = nworker, memo.size.GB = memo.size.GB,
          weighted = weighted, exclude.consp = exclude.consp,
          rand = rand, output.bMNTD = output.bMNTD,
          sig.index=sig.index,
          unit.sum = unit.sum, correct.special = correct.special,
          detail.null=detail.null,
          special.method=special.method,
          ses.cut=ses.cut,rc.cut=rc.cut,conf.cut=conf.cut,
          dirichlet = dirichlet) -> out_icamp


  return(out_icamp)

  detach("package:iCAMP", unload=TRUE)

}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note mntd.obs.z 	 Standardized effect size of MNTD vs. null communities (= (mntd.obs - mntd.rand.mean) / mntd.rand.sd, equivalent to -NTI)
#' @note
#' @note from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/bNTI_Local_Machine.r
#' @return .
#' @export
#' @examples
#' library(phyloseq);library(tidyverse);data("GlobalPatterns")
#' GlobalPatterns %>% subset_samples(SampleType == "Feces") %>% filter_taxa(function(x) sum(x > 0) > 2, TRUE) %>% phyloseq_bNTI_stegen(.,  weighted = TRUE,exclude.conspecifics= FALSE, beta.reps = 4) -> out_phyloseq_bNTI_stegen

phyloseq_bNTI_stegen <- function(physeq, weighted = TRUE,exclude.conspecifics= FALSE, beta.reps = 999){


  ####---------------------- Load R package
  require(picante)

  ####---------------------- Extract data
  # physeq %>%
  #   microbiome::transform(transform = "compositional") -> physeq

  as(otu_table(physeq), "matrix") %>%
    # t() %>%
    as.matrix() -> otu


  physeq@phy_tree %>%
    cophenetic() %>%
    as.matrix() -> dis

  ####----------------------

  ## calculate empirical betaMNTD

  beta.mntd.weighted = as.matrix(picante::comdistnt(t(otu),dis,abundance.weighted=weighted))

  # calculate randomized betaMNTD

  # rand.weighted.bMNTD.comp = array(c(-beta.reps),dim=c(nrow(otu),nrow(otu),beta.reps))
  rand.weighted.bMNTD.comp = array(c(-beta.reps),dim=c(ncol(otu),ncol(otu),beta.reps))
  # dim(rand.weighted.bMNTD.comp)

  for (rep in 1:beta.reps) {

    # rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(otu,taxaShuffle(dis),abundance.weighted=weighted,exclude.conspecifics = exclude.conspecifics))
    rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(otu),taxaShuffle(dis),abundance.weighted=weighted,exclude.conspecifics = exclude.conspecifics))
    print(c(date(),rep))

  }

  weighted.bNTI = matrix(c(NA),nrow=ncol(otu),ncol=ncol(otu))


  for (columns in 1:(ncol(otu)-1)) {
    for (rows in (columns+1):ncol(otu)) {

      rand.vals = rand.weighted.bMNTD.comp[rows,columns,]
      weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals)
      rm("rand.vals")

    }
  }


  rownames(weighted.bNTI) = colnames(otu)
  colnames(weighted.bNTI) = colnames(otu)
  # weighted.bNTI;
  # write.csv(weighted.bNTI,"weighted_bNTI.csv",quote=F);

  # pdf("weighted_bNTI_Histogram.pdf")
  # hist(weighted.bNTI)
  # dev.off()

  return(weighted.bNTI)
  # return(out = list("bNTI" = weighted.bNTI,
  #                   "plot" = hist(weighted.bNTI))
  # )
}



#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note from https://docs.ropensci.org/phylocomr/reference/ph_comdist.html
#' @note species id (same as in the phylogeny, must begin with a letter, not a number or symbol)
#' @note sample_id no spaces, must begin with a letter, not a number or symbol
#' @return .
#' @export
#' @examples
#' library(phyloseq);library(tidyverse);data("GlobalPatterns")
#' GlobalPatterns %>% subset_samples(SampleType == "Feces") %>% filter_taxa(function(x) sum(x > 0) > 3, TRUE) -> physeq; taxa_names(physeq) <- paste0("OTU_", taxa_names(physeq) ); physeq %>% phyloseq_bNTI_phylocomr()
#'

phyloseq_bNTI_phylocomr <- function(physeq,
                                    rand_test = FALSE,
                                    null_model = 0,
                                    randomizations = 999,
                                    abundance = TRUE){

  ####---------------------- Load R package
  require("phylocomr")

  ####---------------------- Extract data
  # physeq %>%
  #   microbiome::transform(transform = "compositional") -> physeq

  as(otu_table(physeq), "matrix") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('sample_id_temp') %>%
    tidyr::pivot_longer(cols = taxa_names(physeq)) %>%
    select(sample_id_temp, value, name) %>%
    mutate(name = as.factor(name)) %>%
    as.data.frame() -> otu


  physeq@phy_tree -> phylo

  ####----------------------

  ph_comdistnt(
    otu,
    phylo,
    rand_test = rand_test,
    null_model = null_model,
    randomizations = randomizations,
    abundance = abundance
  ) -> out

  return(out)

  detach("package:phylocomr", unload=TRUE)

}

#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note
#' @note from https://github.com/FranzKrah/raup_crick/blob/master/raup_crick_abu_par.r
#' @return .
#' @export
#' @examples
#'
#'
#'

phyloseq_raup_crick_abu_par <- function(phyloseq, reps, ncore, classic_metric=FALSE, split_ties=TRUE){



  ####---------------------- Load R package
  require("parallel");  require("doSNOW")

  ####---------------------- Extract data
  # physeq %>%
  #   microbiome::transform(transform = "compositional") -> physeq

  as(otu_table(physeq), "matrix") %>%
    t() %>%
    as.matrix() -> com

  ####----------------------


  pb <- txtProgressBar(max =reps, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  cl <- makeCluster(ncore)
  registerDoSNOW(cl)

  bray.rand <- foreach(randomize = 1:reps,
                       .options.snow = opts,
                       .packages = c("vegan", "picante")) %dopar% {


                         null.dist <- com*0

                         for(i in 1:nrow(com)){

                           com.pa <- (com>0)*1
                           gamma<-ncol(com)
                           occur<-apply(com>0, MARGIN=2, FUN=sum)
                           abundance<-apply(com, MARGIN=2, FUN=sum)
                           com1 <- rep(0,gamma)

                           com1[sample(1:gamma, sum(com.pa[i,]), replace=FALSE, prob=occur)]<-1
                           com1.samp.sp = sample(which(com1>0), (sum(com[i,])-sum(com1)),
                                                 replace=TRUE,prob=abundance[which(com1>0)]);
                           com1.samp.sp = cbind(com1.samp.sp,1)
                           com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum))
                           colnames(com1.sp.counts) = 'counts'
                           com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts))
                           com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts
                           x <- com1
                           null.dist[i,] <- x
                           rm('com1.samp.sp','com1.sp.counts')
                         }
                         as.matrix(vegdist(null.dist, "bray"))
                       }
  stopCluster(cl)

  ## Calculate beta-diversity for obs metacommunity
  bray.obs <- as.matrix(vegdist(com, "bray"))

  ##how many null observations is the observed value tied with?
  null_bray_curtis <- bray.rand
  num_exact_matching_in_null <- lapply(null_bray_curtis, function(x) x==bray.obs)
  num_exact_matching_in_null <- apply(simplify2array(num_exact_matching_in_null), 1:2, sum)

  ##how many null values are smaller than the observed *dissimilarity*?
  num_less_than_in_null <- lapply(null_bray_curtis, function(x) (x<bray.obs)*1)
  num_less_than_in_null <- apply(simplify2array(num_less_than_in_null), 1:2, sum)


  rc = (num_less_than_in_null)/reps; # rc;

  if(split_ties){

    rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
  };


  if(!classic_metric){

    ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1

    rc = (rc-.5)*2
  };

  return(rc)

}



#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note so far sample wise but can be fine tuned using: phyloseq::merge_samples()
#' @note parrallel version: https://stackoverflow.com/questions/46532657/convert-r-apply-statement-to-lapply-for-parallel-processing
#' @note see: https://www.biorxiv.org/content/10.1101/2020.10.15.341966v1.full.pdf
#' @note TODO: # Now generate an object that combines the replicates at each level and time point
#' @return
#' @export
#' @examples
#' require(tidyverse); require(phyloseq); data("GlobalPatterns"); GlobalPatterns %>% phyloseq::merge_samples("SampleType", fun = sum) %>%   phyloseq_kraft_null_model(nsim = 2, meta_sel = sample_variable(.), filtering_expr = "is.na('sample_A')") -> out
#'


phyloseq_kraft_null_model <- function(phyloseq = GlobalPatterns,
                                      meta_sel = c("SampleType",sample_variables(phyloseq)),
                                      nsim = 999,
                                      filtering_expr = "SampleType_A == 'Freshwater (creek)' | SampleType_B == 'Sediment (estuary)'",
                                      verbose = TRUE)
{
  ####---------------------- Load R package -------

  require(phyloseq); require(tidyverse)#; require()

  ####---------------------- Extract data -------

  phyloseq %>%
    filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> phyloseq

  as(otu_table(phyloseq), "matrix") %>%
    data.frame() %>%
    rownames_to_column('Taxa') %>%
    pivot_longer(sample_names(phyloseq)) %>%
    filter(value > 0) -> df

  ####---------------------- old way -------

  # df <- tibble(spp=rep(df$Taxa,df$value), transect=rep(df$name,df$value))

  ####---------------------- XXXXXXXXXXX -------

  phyloseq %>%
    sample_names() %>%
    as.vector() -> samples_names

  ####---------------------- XXXXXXXXXXX -------

  tibble::tibble(Sample_A = samples_names,
                 Sample_B = samples_names) %>%
    tidyr::expand(Sample_A, Sample_B) %>%
    dplyr::filter(Sample_A != Sample_B) -> sample_pw  # remove self comparaisons

  #### ------------- Extracting metadata based on meta_sel parameter

  phyloseq %>%
    sample_data() %>%
    data.frame() %>%
    dplyr::select(one_of(meta_sel)) %>%
    rownames_to_column("sample_id_tmp") -> meta_data

  meta_data_A <- meta_data
  meta_data_B <- meta_data

  names(meta_data_A) <- paste0(names(meta_data), "_A")
  names(meta_data_B) <- paste0(names(meta_data), "_B")

  #### ------------- Joining metadata with pairwise sample comparaisons

  sample_pw %>%
    left_join(meta_data_A,
              by = c("Sample_A" = "sample_id_tmp_A")) %>%
    left_join(meta_data_B,
              by = c("Sample_B" = "sample_id_tmp_B")) %>%
    # filter("SampleType_B" == "Mock")
    filter(rlang::eval_tidy(rlang::parse_expr(filtering_expr))) -> sample_pw_meta


  # Stab_Treat_A  %in% "STAB" & # Sample_A has to be STAB  (Stab_Treat_A colum)
  #        Stab_Treat_B == "TREAT" & # and Sample_B has to be TREAT (Stab_Treat_B colum)
  #        Period_A == Period_B, # and Period_A (of Sample_A has to be the same as the period of Sample_B - no between period ratios)
  #      Treatment_A == Treatment_B) -> sample_pw_meta # only within treatment ratios


  ####---------------------- XXXXXXXXXXX -------

  out_list <- list()

  for (i in 1:nrow(sample_pw_meta)){

    # print(i)

    # df %>%
    #   filter(transect == sample_pw$Sample_A[i] | transect == sample_pw$Sample_B[i]) -> out_list[[i]]
    #
    # names(out_list)[[i]] <- paste0(sample_pw$Sample_A[i], "_", sample_pw$Sample_B[i])

    df %>%
      filter( name %in% c(sample_pw$Sample_A[i],
                          sample_pw$Sample_B[i])) %>%
      filter(value > 0 ) -> out_tmp


    tibble(spp=rep(out_tmp$Taxa,out_tmp$value), transect=rep(out_tmp$name,out_tmp$value)) -> out_list[[i]]

    names(out_list)[[i]] <- paste0(sample_pw$Sample_A[i], "_", sample_pw$Sample_B[i])

  }


  ####---------------------- XXXXXXXXXXX -------


  # out_list[[1]] %>%
  #   ses.beta.function(Nsim = 2)
  #
  # mini_list <- list(  out_list[[1]],   out_list[[2]],   out_list[[3]])

  # names(final_out) <- names(out_list)
  #
  # mini_list %>%
  #   lapply(ses.beta.function) -> final_out

  final_out <- list()

  out_list %>%
    lapply(ses.beta.function) %>%
    bind_rows(.) %>%
    mutate(gamma = as.double(gamma),
           obs.mean.alpha = as.double(obs.mean.alpha),
           obs.beta  = as.double(obs.beta ),
           mean.null.beta = as.double(mean.null.beta),
           sd.null.beta  = as.double(sd.null.beta ),
           ses.beta   = as.double(ses.beta)) %>%
    mutate(SI = 1 - ((abs(obs.beta - mean.null.beta))/obs.beta )) %>%
    mutate(DS = (abs(obs.beta - mean.null.beta))/obs.beta) -> final_out

  ####---------------------- XXXXXXXXXXX -------

  ses.beta.function <- function(gdata, Nsim=nsim)
    # 'gdata' is the    data file, 'Nsim' is the number of randomizations for calculating expected/null beta
  {
    ####---------------------- XXXXXXXXXXX -------

    require(tidyverse)

    ####---------------------- XXXXXXXXXXX -------

    gdata %>%
      mutate(transect = as.vector(transect)) %>%
      mutate(spp = as.vector(gdata$spp)) -> gdata

    ####---------------------- XXXXXXXXXXX -------

    sample_A = unique(gdata$transect)[1]
    sample_B = unique(gdata$transect)[2]

    if(verbose == TRUE){
      print(paste0("processing ", sample_A, " and ", sample_B ))
    }

    ####---------------------- XXXXXXXXXXX -------

    plot.gamma=length(unique(gdata$spp)) #calculate the total number of species at the site
    transect.spp=tapply(gdata$spp,gdata$transect,unique) #generate species list for each transect
    obs.mean.alpha=mean(sapply(transect.spp,length)) #calculate average number of species per transect
    obs.beta=1-obs.mean.alpha/plot.gamma #calculate observed beta partition
    rand.mean.alpha=vector(length=Nsim) #create empty vector to be filled in with randomly generated alpha values
    for(j in 1:Nsim) #start loop for simulations
    {
      if(verbose == TRUE){

        print(paste0("simulation:  ", j ))
      }

      samp=sample(gdata$transect,length(gdata$transect),replace=F)

      #swaps order of plotnames
      swap.data=data.frame("transect"=samp,"spp"=gdata$spp)

      #assigns random plotnames to individuals
      rand.transect.spp=tapply(swap.data$spp,swap.data$transect,unique)

      #generate species list for each transect
      rand.mean.alpha[j]=mean(sapply(rand.transect.spp,length))

      #calculate average number of species per transect
    } #end loop for simulations
    null.plot.beta=1-rand.mean.alpha/plot.gamma #calculates the 1000 random beta values

    mean.null.beta=mean(null.plot.beta) #calculates the mean of the random beta values

    sd.null.beta=sd(null.plot.beta) #calculates the sd of the random beta values

    ses.beta=(obs.beta-mean.null.beta)/sd.null.beta #calculates the deviation of the observed from expected (random) beta

    out <- c("sample_A" =sample_A, "sample_B" = sample_B,"gamma"=plot.gamma,"obs.mean.alpha"=obs.mean.alpha,"obs.beta"=obs.beta,"mean.null.beta"=mean.null.beta,"sd.null.beta"=sd.null.beta,"ses.beta"=ses.beta)

    # # rownames(out)
    # unique(gdata$transect) %>%  as.character() -> samples
    #
    # paste0(samples[1], "_vs_" ,samples[2])
    #
    # lapply(strsplit(unique(gdata$transect, ","), as.character)
    #
    #        -> opt$raw_file_pattern

    return(out)
  }

  ####---------------------- XXXXXXXXXXX -------

  return(final_out)
}

