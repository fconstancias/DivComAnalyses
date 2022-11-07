#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note see: https://chiliubio.github.io/microeco_tutorial/model-based-class.html#trans_nullmodel-class
#' @note
#' @note .
#' @return .
#' @export
#' @examples
#'
#'require(phyloseq); require(tidyverse)
#'data("GlobalPatterns")
#'GlobalPatterns -> physeq
#'
#'set.seed(12344566)
#'sample_data(physeq)$pH <- rnorm(nsamples(physeq), mean=6, sd=2)
#'set.seed(12344566)
#'sample_data(physeq)$temp  <- rnorm(nsamples(physeq), mean=22, sd=6)
#'physeq %>% subset_samples(SampleType == "Feces") %>% filter_taxa(function(x) sum(x > 0) > 2, TRUE) -> physeq
#'
#'physeq %>% phyloseq_null_model_microeco(., env_cols = c("pH", "temp"),  method = c("mantel_corr", "betaMNTD",  "RCbray",  "cal_NTI", "cal_process"), test = TRUE) -> out
#'

readRDS((url("https://github.com/fconstancias/DivComAnalyses/blob/9e6ebc69e0c5e136c910186a7136e51af4cecc13/data-raw/ps_invivo.RDS?raw=true" ))) %>% filter_taxa(function(x) sum(x > 0) > 1, TRUE) %>% phyloseq_null_model_microeco(env_cols = NULL, c("mantel_corr", "betaMNTD", "cal_NTI","RCbray","cal_process")) -> out

phyloseq_null_model_microeco <- function(physeq,
                                         method = c("mantel_corr", "betaMPD", "betaMNTD", "RCbray", "cal_NRI", "cal_NTI", "cal_Cscore", "cal_tNST", "cal_process"),
                                         env_cols = c(NULL, "pH", "temp"),
                                         mantel_r.type = "pearson",
                                         mantel_nperm = 999,
                                         mantel_mult = "holm",
                                         abundance_weighted = TRUE,
                                         exclude_conspecifics = FALSE,
                                         runs = 1000,
                                         iterations = 1000,
                                         null.model = c("taxa.labels"),
                                         Cscore_group = NULL,
                                         tNST_group = "Group",
                                         Cscore_normalise = FALSE,
                                         tNST_group_dist = "bray",
                                         threads = 4,
                                         test = FALSE)
{

  ####---------------------- Load R package

  require(microeco); require(phyloseq); require(file2meco); require(tidyverse)

  ####---------------------- Extract data

  # data(dataset)
  # # data(env_data_16S)
  #
  # dataset %>%
  #   file2meco::meco2phyloseq() -> physeq
  #
  #
  # data("GlobalPatterns")
  # # GlobalPatterns -> physeq
  #
  # set.seed(12344566)
  # sample_data(physeq)$pH <- rnorm(nsamples(physeq), mean=6, sd=2)
  # set.seed(12344566)
  # sample_data(physeq)$temp  <- rnorm(nsamples(physeq), mean=22, sd=6)

  physeq %>%
    file2meco::phyloseq2meco(.) -> data

  data$sample_table -> env_data

  ####---------------------- null model analysis - trans_nullmodel {microeco}

  t1 <- trans_nullmodel$new(data,
                            taxa_number = ifelse(test == TRUE, 20, ntaxa(physeq)),
                            # env_cols = env_cols,
                            add_data = env_data)#ifelse(env_cols == NULL, "NULL", env_data) )

  out = NULL

  ####---------------------- mantel_corr
  if ("mantel_corr"  %in% method)
  {
    t1$cal_mantel_corr(use_env = env_cols,
                       break.pts = seq(0, 1, 0.02),
                       r.type = mantel_r.type,
                       nperm = mantel_nperm,
                       mult = mantel_mult,
                       progressive = TRUE)

    t1$plot_mantel_corr() -> mantel_cor

    out$mantel_cor <- mantel_cor
  }
  ####---------------------- betaMPD
  #Calculate betaMPD (mean pairwise distance). Same with picante::comdist function, but faster.

  if ("betaMPD"  %in% method)
  {
    t1$cal_betampd(abundance.weighted = abundance_weighted)

    t1$res_betampd %>%
      data.frame() -> res_betampd

    t1$cal_ses_betampd(
      runs = ifelse(test == TRUE, 9, runs) ,
      null.model = null.model,
      abundance.weighted = abundance_weighted,
      iterations = ifelse(test == TRUE, 9, iterations))

    t1$res_ses_betampd %>%
      data.frame() -> res_ses_betampd

    out$res_betampd <- res_betampd
    out$res_ses_betampd <- res_ses_betampd
  }

  ####---------------------- betaMNTD
  #Calculate betaMNTD (mean nearest taxon distance). Same with picante::comdistnt package, but faster.

  if ("betaMNTD"  %in% method)
  {
    t1$cal_betamntd(abundance.weighted = abundance_weighted,
                    use_iCAMP = FALSE, use_iCAMP_force = TRUE)

    t1$res_betamntd %>%
      data.frame() -> res_betamntd

    t1$cal_ses_betamntd(
      runs = ifelse(test == TRUE, 9, runs) ,
      null.model = null.model,
      nworker = threads,
      exclude.conspecifics = exclude_conspecifics,
      abundance.weighted = abundance_weighted,
      iterations = ifelse(test == TRUE, 9, iterations))

    t1$res_ses_betamntd %>%
      data.frame() -> res_ses_betamntd

    out$res_betamntd <- res_betamntd
    out$res_ses_betamntd <- res_ses_betamntd
  }

  ####---------------------- Raup–Crick
  #Calculate Bray–Curtis-based Raup–Crick (RCbray).

  if ("RCbray"  %in%  method)
  {
    t1$cal_rcbray(
      runs = ifelse(test == TRUE, 2, runs) ,
      verbose = TRUE,
      null.model = "independentswap")

    t1$res_rcbray ->  out$res_rcbray
  }

  ####---------------------- cal_NRI
  # Calculates Nearest Relative Index (NRI), equivalent to -1 times the standardized effect size of MPD.
  if ("cal_NRI"  %in%  method)
  {
    t1$cal_NRI(
      null.model = null.model,
      abundance.weighted = abundance_weighted,
      runs = ifelse(test == TRUE, 2, runs))

    t1$res_NRI ->  out$res_NRI
  }

  ####---------------------- cal_NTI
  # Calculates Nearest Taxon Index (NTI), equivalent to -1 times the standardized effect size of MNTD.

  if ("cal_NTI"  %in%  method)
  {
    t1$cal_NTI(
      null.model = null.model,
      abundance.weighted = abundance_weighted,
      runs = ifelse(test == TRUE, 2, runs))

    t1$res_NTI -> out$res_NTI
  }

  ####---------------------- cal_Cscore
  # Calculates the (normalised) mean number of checkerboard combinations (C-score) using C.score function in bipartite package.

  if ("cal_Cscore"  %in%  method)
  {
    t1$cal_Cscore(by_group = Cscore_group,
                  normalise = Cscore_normalise) -> Cscore

    out$Cscore <- Cscore
  }

  ####---------------------- cal_tNST
  # Calculate normalized stochasticity ratio (NST) based on the tNST function of NST package.

  if ("cal_tNST"  %in%  method)
  {
    t1$cal_tNST(group = tNST_group,
                dist.method = tNST_dist,
                output.rand = TRUE,
                SES = TRUE)

    t1$res_tNST -> out$res_tNST

    t1$cal_tNST_test(method = "nst.boot")
  }

  ####---------------------- cal_process
  # Infer the ecological processes according to ses.betaMNTD ses.betaMPD and rcbray.

  if ("cal_process"  %in%  method)
  {
    t1$cal_process(use_betamntd = TRUE)

    t1$res_process -> out$res_process
  }

  ####---------------------- return output
  return(out)
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note see: https://chiliubio.github.io/microeco_tutorial/explainable-class.html#trans_env-class
#' @note
#' @note .
#' @return .
#' @export
#' @examples
#'
#'require(phyloseq); require(tidyverse)
#'readRDS((url("https://github.com/fconstancias/DivComAnalyses/blob/9e6ebc69e0c5e136c910186a7136e51af4cecc13/data-raw/ps_invivo.RDS?raw=true" ))) %>% filter_taxa(function(x) sum(x > 0) > 110, TRUE) %>% phyloseq_env_microeco(env_cols = c("BW", "merged_pc"), cal_diff_group = "treatment_grouped")

phyloseq_env_microeco <- function(physeq, method = c("cal_diff", "cal_autocor"), env_cols = NULL, cal_diff_group = "Group", cal_diff_m = "wilcox", cal_ordination_m = c("dbRDA", "CCA", "RDA"), ordination_feature_sel = FALSE){


  ####---------------------- Load R package

  require(microeco); require(phyloseq); require(file2meco); require(tidyverse)

  ####---------------------- Extract data

  # data("GlobalPatterns")
  # GlobalPatterns  %>%  subset_samples(SampleType %in% c("Feces",  "Soil")) -> physeq
  #
  # set.seed(12344566)
  # sample_data(physeq)$pH <- rnorm(nsamples(physeq), mean=6, sd=2)
  # set.seed(12344566)
  # sample_data(physeq)$temp  <- rnorm(nsamples(physeq), mean=22, sd=6)

  physeq %>%
    file2meco::phyloseq2meco(.) -> data

  data$sample_table -> env_data

  ####---------------------- null model analysis - trans_nullmodel {microeco}

  t1 <- trans_env$new(dataset = data,
                      env_cols = env_cols,
                      add_data = NULL,
                      character2numeric = TRUE,
                      complete_na = FALSE)

  out = NULL


  ####---------------------- cal_diff

  if ("cal_diff"  %in% method)
  {

    t1$cal_diff(group = cal_diff_group, method = cal_diff_m, anova_set = NULL,  by_group = NULL, p_adjust_method = "fdr")

    t1$res_diff -> out$res_diff

    res_diff_plot <- vector("list", length(env_cols))
    names(res_diff_plot) <- env_cols

    tmp <- list()
    for(i in colnames(t1$data_env)){
      tmp[[i]] <- t1$plot_diff(measure = i, add_sig_text_size = 5, xtext_size = 12) + theme(plot.margin = unit(c(0.1, 0, 0, 1), "cm"))
    }

    # for (i in names((res_diff_plot))
    # {
    # t1$plot_diff(color_values = RColorBrewer::brewer.pal(env_data[,cal_diff_group] %>%  unique() %>%  length(), "Dark2"),
    #              measure = env_cols[i],
    #              group = cal_diff_group,
    #              add_sig = TRUE,
    #              add_sig_label = "Significance",
    #              add_sig_text_size = 3.88,
    #              use_boxplot = TRUE,
    #              boxplot_add = "jitter",
    #              order_x_mean = FALSE,
    #              y_start = 1.01,
    #              y_increase = 0.05,
    #              xtext_angle = NULL,
    #              xtext_size = 15,
    #              ytitle_size = 17,
    #              barwidth = 0.9) -> res_diff_plot[[i]]
    # }

    out$res_diff_plot <- tmp
  }

  ####---------------------- cal_autocor

  if ("cal_autocor"  %in% method)
  {

    t1$cal_autocor(
      group = cal_diff_group,
      color_values = RColorBrewer::brewer.pal(env_data[,cal_diff_group] %>%  unique() %>%  length(), "Dark2"),
      alpha = 0.8,
      upper = list(continuous = GGally::wrap("cor", method= "spearman"))) -> out$cal_autocor

  }

  ####---------------------- cal_autocor


  if ("cal_ordination"  %in% method)
  {

    t1$cal_ordination(
      method = cal_ordination_m,
      feature_sel = ordination_feature_sel,
      taxa_level = NULL,
      taxa_filter_thres = NULL,
      use_measure = NULL,
      add_matrix = NULL)

  }


  ####---------------------- cal_mantel


  if ("cal_mantel"  %in% method)
  {

    t1$cal_mantel(
      select_env_data = env_cols,
      partial_mantel = TRUE,
      add_matrix = NULL,
      use_measure = "bray",
      method = "pearson",
      p_adjust_method = "fdr")

  }

  return(out)
}




#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note see: https://chiliubio.github.io/microeco_tutorial/explainable-class.html#trans_func-class
#' @note
#' @note .
#' @return .
#' @export
#' @examples
#'
#'require(phyloseq); require(tidyverse)
#'data("GlobalPatterns")
#'GlobalPatterns -> physeq
#'
#'set.seed(12344566)
#'sample_data(physeq)$pH <- rnorm(nsamples(physeq), mean=6, sd=2)
#'set.seed(12344566)
#'sample_data(physeq)$temp  <- rnorm(nsamples(physeq), mean=22, sd=6)
#'physeq %>% subset_samples(SampleType == "Feces") %>% filter_taxa(function(x) sum(x > 0) > 2, TRUE) -> physeq
#'
#'physeq %>% phyloseq_null_model_microeco(., env_cols = c("pH", "temp"),  method = c("mantel_corr", "betaMNTD",  "RCbray",  "cal_NTI", "cal_process"), test = TRUE) -> out
#'
#'
readRDS((url("https://github.com/fconstancias/DivComAnalyses/blob/9e6ebc69e0c5e136c910186a7136e51af4cecc13/data-raw/ps_invivo.RDS?raw=true" ))) %>% filter_taxa(function(x) sum(x > 0) > 110, TRUE) %>% phyloseq_func_microeco(prok_database = c("FAPROTAX")) -> out

phyloseq_func_microeco <- function(physeq,method = c("cal_spe_func", "cal_tax4fun2"), prok_database = c("FAPROTAX", "NJC19"), fungi_database = c("FUNGuild", "FungalTraits"), abundance_weighted = TRUE, blast_tool_path = "~/Documents/ncbi-blast-2.13.0+/bin/", path_to_reference_data = "~/Documents/Ref99NR/Tax4Fun2_ReferenceData_v2/", num_threads = 4){


  ####---------------------- Load R package

  require(microeco); require(phyloseq); require(file2meco); require(tidyverse)

  ####---------------------- Extract data

  # data("GlobalPatterns")
  # GlobalPatterns  %>%  subset_samples(SampleType %in% c("Feces",  "Soil")) -> physeq
  #
  # set.seed(12344566)
  # sample_data(physeq)$pH <- rnorm(nsamples(physeq), mean=6, sd=2)
  # set.seed(12344566)
  # sample_data(physeq)$temp  <- rnorm(nsamples(physeq), mean=22, sd=6)

  # data("dataset")
  #
  # dataset %>%
  #   file2meco::meco2phyloseq() -> pss

  # tax_table(physeq) <-  gsub(tax_table(physeq), pattern = "unknown", replacement = "s__")

  otu_table_trans <- as.data.frame(physeq@otu_table@.Data, check.names = FALSE, stringsAsFactors = FALSE)
  sample_table_trans <- data.frame(phyloseq::sample_data(physeq), check.names = FALSE, stringsAsFactors = FALSE)
  tax_table_trans <- as.data.frame(physeq@tax_table@.Data, check.names = FALSE, stringsAsFactors = FALSE)
  tax_table_trans %>%  tidy_taxonomy() -> tax_table_trans
  phylo_tree_trans <- physeq@phy_tree
  seq_trans <- physeq@refseq

  data <- microtable$new(sample_table = sample_table_trans, otu_table = otu_table_trans,
                            tax_table = tax_table_trans, phylo_tree = phylo_tree_trans, rep_fasta = seq_trans)

  # physeq %>%
    # file2meco::phyloseq2meco(physeq = .) -> data

  # data@rep_fasta <- physeq@refseq

  data$sample_table -> env_data

  ####---------------------- null model analysis - trans_nullmodel {microeco}

  t1 <- trans_func$new(data)

  out = NULL

  ####---------------------- cal_diff

  if ("cal_spe_func"  %in% method)
  {

    t1$cal_spe_func(prok_database = prok_database)

    t1$res_spe_func -> out$res_spe_func

    t1$cal_spe_func_perc(abundance_weighted = abundance_weighted)

    t1$res_spe_func_perc -> out$res_spe_func_perc

    t1$show_prok_func() -> tt

    t1$plot_spe_func_perc(
      filter_func = NULL,
      use_group_list = TRUE,
      add_facet = TRUE,
      order_x = NULL,
      color_gradient_low = "#00008B",
      color_gradient_high = "#9E0142"
    ) -> out$plot_spe_func_perc

    # t1$print

  }

  ####---------------------- cal_diff

  if ("cal_tax4fun2"  %in% method)
  {

    t1$cal_tax4fun2(blast_tool_path = blast_tool_path,
                    path_to_reference_data = path_to_reference_data,
                    normalize_by_copy_number = T,
                    min_identity_to_reference = 97,
                    use_uproc = T,
                    num_threads = num_threads,
                    normalize_pathways = F)

    tax4fun2 <- microtable$new(otu_table = t1$res_tax4fun2_pathway, tax_table = Tax4Fun2_KEGG$ptw_desc, sample_table = dataset$sample_table)

    tax4fun2$tidy_dataset()

    tax4fun2$cal_abund()

    func2 <- trans_abund$new(tax4fun2, taxrank = "Level.2")#, groupmean = "treatment_grouped")
    func2$plot_bar(legend_text_italic = FALSE) -> fun_plot

    # calculate functional redundancies
    t1$cal_tax4fun2_FRI()

    t1$res_tax4fun2_aFRI -> out$res_tax4fun2_aFRI

    t1$res_tax4fun2_rFRI -> out$res_tax4fun2_rFRI

    # t1$res_tax4fun2_otu_table_reduced_aggregated

    # t1$cal_abund()

    # t1$plot_spe_func_perc(
    #   filter_func = NULL,
    #   use_group_list = TRUE,
    #   add_facet = TRUE,
    #   order_x = NULL,
    #   color_gradient_low = "#00008B",
    #   color_gradient_high = "#9E0142"
    # )
  }

  return(out)
}



