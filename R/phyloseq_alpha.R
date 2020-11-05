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
#' library(phyloseq)
#' library(tidyverse)
#' data("GlobalPatterns")
#'
#' phyloseq_alphas(GlobalPatterns, phylo = TRUE) -> results
#'

phyloseq_alphas <- function(physeq,
                            phylo)
{
  # options(getClass.msg=FALSE)
  require(metagMisc)
  require(microbiome)
  if (phylo == TRUE)
  {

    bind_cols(
      microbiome::meta(physeq) %>%
        rownames_to_column('id1'),
      microbiome::alpha(physeq, index = c("diversity_inverse_simpson",
                                          "diversity_shannon",
                                          "evenness_pielou",
                                          "diversity_coverage")) %>%
        rownames_to_column('id2'),
      phyloseq::estimate_richness(physeq,measures =c("Observed", "Chao1", "ACE")) %>%
        as.data.frame() %>%
        rownames_to_column('id3'),
      metagMisc::phyloseq_phylo_div(physeq) %>%
        rownames_to_column('id4')) -> tmp


    metagMisc::phyloseq_phylo_ses(physeq,
                                  package = "PhyloMeasures",
                                  verbose = FALSE) -> tmp2

    left_join(tmp,
              tmp2,
              by = c("id1" = "SampleID"),
              suffix = c("", "_PM"),
              keep = TRUE ) -> tmp3

tmp3 %>%
  mutate(bNTI = (MNTD_PM - MNTD.rand.mean) / MNTD.rand.sd ) -> alpha

return(alpha)
  }else
  {

    bind_cols(
      microbiome::meta(physeq) %>%
        rownames_to_column('id'),
      microbiome::alpha(physeq) %>%
        rownames_to_column('id'),
      phyloseq::estimate_richness(physeq,measures =c("Observed", "Chao1", "ACE")) %>%
        as.data.frame() %>%
        rownames_to_column('id')
    ) -> alpha

    return(alpha)
    detach("package:metagMisc", unload=TRUE);     detach("package:microbiome", unload=TRUE)

  }

}

#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .remove ASV assigned to Chloroplast and Mitochondria
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'


plot_alphas <- function(alphas,
                        measure = c("Observed", "diversity_shannon"),
                        x_group,
                        colour_group,
                        fill_group,
                        shape_group,
                        facet_group = FALSE,
                        test_group,
                        test_group_2){
  alphas %>%
    pivot_longer(cols = all_of(measure), values_to = "value", names_to = 'alphadiversiy', values_drop_na  = TRUE) %>%
    mutate(alphadiversiy = fct_relevel(alphadiversiy, measure)) %>%
    ggplot(aes_string(x_group,
                      "value",
                      colour = colour_group,
                      fill = fill_group)) +
    geom_boxplot(outlier.colour = NA, alpha=0.7) +
    # ggbeeswarm::geom_beeswarm(size=1, alpha=0.2,
    #                           position=position_jitterdodge(dodge.width=0.9)) +
    geom_jitter(size=1, position = position_jitterdodge(dodge.width=1),
                aes_string(shape = shape_group)) +
    # geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
    ylab("Diversity index")  + xlab(NULL) + theme_light() -> p

  p + facet_grid(as.formula(paste0("alphadiversiy ~ ",paste(facet_group))), scales = "free_y", space = "fixed") -> p 
  
ggpubr::compare_means(formula = as.formula(paste0("value ~ ", paste0(test_group))),
                        group.by = c("alphadiversiy", test_group_2),
                        data = p$data,
                        method = "wilcox.test",
                        p.adjust.method = "fdr") %>%
    select(-.y., -p.format, -p.signif) %>%
    arrange(p) %>%
    mutate(signif = ifelse(p.adj <= 0.05, 'SIGN', 'NS')) -> stat


  out <- list("plot" = p,
              "stat" = stat)

  return(out)

}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .remove ASV assigned to Chloroplast and Mitochondria
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
correlate_alpha <- function(alphas_all,
                            colnums_to_plot,
                            colour,
                            shape,
                            method)
{
  # define a function to plot scatter plot
  my_fn <- function(data, mapping, ...){
    p <- ggplot(data = data, mapping = mapping) +
      geom_point() +
      geom_smooth(method=lm, ...)
    p
  }

  # where alphas_all is a dataframe with alpha div and covaraibles i.e., Observed", "evenness_pielou",
  # "diversity_shannon", "SES.MPD","TMAO", "Choline" and also here diet as columns

  alphas_all %>%
    GGally::ggpairs(columns = colnums_to_plot,
                   # ggplot2::aes(colour = colour, shape = shape ),
                   # legend = 1,
                    progress = FALSE,
                    upper = list(
                      continuous = GGally::wrap('cor', method = method)
                    ),
                    lower = list(continuous = my_fn)) -> corplot

  return(corplot)
}
