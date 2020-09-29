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


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note . correlates ASV/Taxa with metadata
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'library(phyloseq)
#'data("enterotype")
#'
#'
#'
#'
#' #"spearman"
#'
phyloseq_correlate_taxa <- function(ps_tmp,
                                    grouping_column,
                                    adjustment= 3,
                                    cor_variables,
                                    method)
{
  ps_tmp %>%
    tax_table() %>%
    data.frame() -> tmp

  ps_tmp %>%
    transform_sample_counts(function(x) x/sum(x) * 100) %>%
    microbiome::transform("log10") %>%
    microbiomeSeq::taxa.env.correlation(grouping_column="run", method= method, pvalue.threshold=0.05,
                                        padjust.method="fdr", adjustment=3, num.taxa=20, select.variables = cor_variables) -> env.taxa.cor

  # plot
  p <- microbiomeSeq::plot_taxa_env(env.taxa.cor)

  p$data %>%
    dplyr::left_join(tmp %>% rownames_to_column("ASV"),
                     by = c("Taxa" = "ASV")) %>%
    dplyr::select(-Taxa) %>%
    dplyr::rename(Taxa = Strain) %>%
    microbiomeSeq::plot_taxa_env() +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab") -> plot

  return(list("plot" = plot ,
              "table" = env.taxa.cor))
}
