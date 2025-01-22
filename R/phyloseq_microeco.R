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

# readRDS((url("https://github.com/fconstancias/DivComAnalyses/blob/9e6ebc69e0c5e136c910186a7136e51af4cecc13/data-raw/ps_invivo.RDS?raw=true" ))) %>% filter_taxa(function(x) sum(x > 0) > 1, TRUE) %>% phyloseq_null_model_microeco(env_cols = NULL, c("mantel_corr", "betaMNTD", "cal_NTI","RCbray","cal_process")) -> out

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

phyloseq_env <- function(physeq, method = c("cal_diff", "cal_autocor"), env_cols = NULL, cor_method = "spearman", cal_diff_group = "Group", cal_diff_m = "wilcox", cal_ordination_m = c("dbRDA", "CCA", "RDA"), ordination_feature_sel = FALSE, color_values = NULL, use_data = c("Genus", "all", "other"), p_adjust_method = "none", p_adjust_type = "Type"){
  
  
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
      color_values = ifelse(color_values == NULL,
                            RColorBrewer::brewer.pal(env_data[,cal_diff_group] %>%  unique() %>%  length(), "Dark2"),
                            color_values),
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
  
  if ("cal_cor"  %in% method)
  {
    # create trans_env object
    # t1 <- trans_env$new(dataset = dataset, add_data = dataset$sample_table, env_cols =  env_cols)
    # calculate correlations
    # t1$cal_abund
    t1$cal_cor(use_data = use_data,p_adjust_method = p_adjust_method, p_adjust_type = p_adjust_type, select_env_data = env_cols,by_group = "group", cor_method = cor_method)
    # plot the correlation heatmap
    t1$plot_cor() -> out$plot_cor
    
    # t1$plot_cor(filter_feature = c("*", "**", "***", "****"))
    
    t1$res_cor -> out$res_cor
    
    t1$res_cor %>%
      filter(AdjPvalue <= 0.05 & Correlation > 0) -> t1$res_cor
    
    t1$plot_cor() -> out$plot_cor
    
    out$plot_cor + scale_fill_gradient2(limits = c(-1, 1), low = "#2C7BB6",
                                        mid = "white", high = "#D7191C") ->  out$plot_cor
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
# readRDS((url("https://github.com/fconstancias/DivComAnalyses/blob/9e6ebc69e0c5e136c910186a7136e51af4cecc13/data-raw/ps_invivo.RDS?raw=true" ))) %>% filter_taxa(function(x) sum(x > 0) > 110, TRUE) %>% phyloseq_func_microeco(prok_database = c("FAPROTAX")) -> out

phyloseq_func <- function(physeq,method = c("cal_spe_func", "cal_tax4fun2"), prok_database = c("FAPROTAX", "NJC19"), fungi_database = c("FUNGuild", "FungalTraits"), abundance_weighted = TRUE, blast_tool_path = "~/Documents/ncbi-blast-2.13.0+/bin/", path_to_reference_data = "~/Documents/Ref99NR/Tax4Fun2_ReferenceData_v2/",  min_identity_to_reference = 97, num_threads = 4){
  
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
    
    t1$show_prok_func() -> out$show_prok_func
    
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
                    min_identity_to_reference = min_identity_to_reference,
                    use_uproc = T,
                    num_threads = num_threads,
                    normalize_pathways = F)
    
    tax4fun2 <- microtable$new(otu_table = t1$res_tax4fun2_pathway, tax_table = Tax4Fun2_KEGG$ptw_desc, sample_table = env_data)
    
    tax4fun2$tidy_dataset()
    
    tax4fun2$cal_abund() #select_cols = "Level.3")
    
    # tax4fun2$taxa_abund()
    
    func2 <- trans_abund$new(tax4fun2, taxrank = "Level.2")#, groupmean = "treatment_grouped")
    func2$plot_bar(legend_text_italic = FALSE) -> out$fun_plot_l2
    
    func1 <- trans_abund$new(tax4fun2, taxrank = "Level.1")#, groupmean = "treatment_grouped")
    func1$plot_bar(legend_text_italic = FALSE) -> out$fun_plot_l1
    
    func3 <- trans_abund$new(tax4fun2, taxrank = "Level.3")#, groupmean = "treatment_grouped")
    func3$plot_bar(legend_text_italic = FALSE) -> out$fun_plot_l3
    
    if ("cal_tax4fun2_FRI"  %in% method)
    {
      # calculate functional redundancies
      t1$cal_tax4fun2_FRI()
      
      t1$res_tax4fun2_aFRI -> out$res_tax4fun2_aFRI
      
      t1$res_tax4fun2_rFRI -> out$res_tax4fun2_rFRI
    }
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


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note see: https://chiliubio.github.io/microeco_tutorial/model-based-class.html
#' @note .
#' @return .
#' @export
#' @examples
#'


phyloseq_diff <- function(physeq = ps_up %>% subset_samples(Sample == "Plaque"), method = "lefse", group = "Time", fix_formula = "Time", alpha = 0.001, lefse_subgroup = NULL, linda_formula = "~ Sample_Type + Time +  (1|Subject)",
                          taxa_level = "all",  filter_thres =  0.00001, p_adjust_method = "fdr", lda_threshold = 2,
                          plot_pal = RColorBrewer::brewer.pal(8, "Dark2"), group_order = NULL, add_sig = FALSE, keep_prefix = TRUE,
                          plot_type2 = "barerrorbar", add_sig_plot2 = FALSE, errorbar_color_black = TRUE,
                          taxa_level_2 = "Species"){
  
  ####---------------------- Load R package
  
  require(microeco); require(phyloseq); require(file2meco); require(tidyverse)
  
  out=NULL
  
  ####---------------------- Extract data
  # sample_data(physeq)$temp  <- rnorm(nsamples(physeq), mean=22, sd=6)
  
  physeq %>%
    file2meco::phyloseq2meco(.) -> data
  
  data$sample_table -> env_data
  
  ####---------------------- feature - trans_diff {microeco}
  
  
  
  
  if ("ancombc2"  %in% method)
  {
    t1 <- trans_diff$new(dataset = data, method = method, group = group, remove_unknown = TRUE,
                         alpha = alpha, lefse_subgroup = lefse_subgroup,
                         taxa_level = taxa_level, filter_thres = filter_thres,
                         p_adjust_method = p_adjust_method)
    
    out$res_diff_raw <- t1$res_diff_raw
    out$abund <- t1$res_abund
    out$res_diff <- t1$res_diff
    
    t1$res_diff %<>% subset(P.adj <= 0.05 & passed_ss == "TRUE")
    
    
    out$diff_bar  <- t1$plot_diff_bar(keep_full_name = FALSE, 
                                      heatmap_cell =  "P.adj",
                                      heatmap_sig = "Significance",
                                      heatmap_x = "Factors",
                                      heatmap_y = "Taxa",
                                      heatmap_lab_fill = "P.adj")
    
    out$diff_bar  <- out$diff_bar + theme(legend.position = "none")
    
    
    # out$diff_abund  <- t1$plot_diff_abund(plot_type = "ggboxplot", select_taxa = t1$plot_diff_bar_taxa , y_start = 0.05, y_increase = 0.1, color_values = plot_pal)
    # 
    # out$diff_abund  <- out$diff_abund + scale_y_sqrt() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank()) + ylab("Proportion (sqrt)") 
    # out$combined_plot <- out$diff_bar%>% aplot::insert_right(out$diff_abund, width = 1)
    # out$combined_plot 
    # 
    
    
  }
  # see t1$res_diff for the result
  # From v0.8.0, threshold is used for the LDA score selection.
  
  if ("lefse"  %in% method)
  {
    t1 <- trans_diff$new(dataset = data, method = method, group = group, remove_unknown = TRUE,
                         alpha = alpha, lefse_subgroup = lefse_subgroup,
                         taxa_level = taxa_level, filter_thres = filter_thres,
                         p_adjust_method = p_adjust_method)
    
    
    out$abund <- t1$res_abund
    out$res_diff <- t1$res_diff
    
    out$diff_bar <- t1$plot_diff_bar(
      threshold = lda_threshold,
      color_values = plot_pal,
      color_group_map = FALSE,
      use_number = use_number,
      select_group = NULL,
      keep_full_name = FALSE,
      keep_prefix = keep_prefix,
      group_order = group_order,
      group_aggre = TRUE,
      group_two_sep = TRUE,
      coord_flip = TRUE,
      add_sig = add_sig,
      add_sig_increase = 0.1,
      add_sig_text_size = 5,
      xtext_angle = 45,
      xtext_size = 10,
      axis_text_y = 12,
      heatmap_cell = ifelse(method %in% c("ancombc2"), "P.adj" ,"P.unadj"),
      heatmap_sig = "Significance",
      heatmap_x = "Factors",
      heatmap_y = "Taxa",
      heatmap_lab_fill = ifelse(method %in% c("ancombc2"), "P.adj" ,"P value"))
    
    out$diff_abund <- t1$plot_diff_abund(group_order = group_order, select_taxa = t1$plot_diff_bar_taxa, plot_type = plot_type2, 
                                         add_sig = add_sig_plot2, errorbar_addpoint = FALSE, errorbar_color_black = errorbar_color_black, 
                                         color_values = plot_pal)
    
    out$diff_bar  <- out$diff_bar + theme(legend.position = "none")
    out$diff_abund  <- out$diff_abund + scale_y_sqrt() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank()) + ylab("Proportion (sqrt)") 
    out$combined_plot <- out$diff_bar%>% aplot::insert_right(out$diff_abund, width = 0.5)
    out$combined_plot 
  }
  
  if ("linda"  %in% method)
  {
    
    t1 <- trans_diff$new(dataset = data, formula = linda_formula , method = "linda", remove_unknown = TRUE,
                         alpha = alpha,
                         taxa_level = taxa_level, filter_thres = filter_thres,
                         p_adjust_method = p_adjust_method)
    
    t1$res_diff %<>% subset(P.adj <= 0.05) # subset(Significance %in% c("*","**","***"))
    
    out$diff_bar  <- t1$plot_diff_bar(keep_full_name = FALSE, 
                                      heatmap_cell =  "P.adj",
                                      heatmap_sig = "Significance",
                                      heatmap_x = "Factors",
                                      heatmap_y = "Taxa",
                                      heatmap_lab_fill = "P.adj")
    
    out$diff_bar  <- out$diff_bar + theme(legend.position = "none")
    
    out$res_diff <- t1$res_diff
    # y_start and y_increase control the position of labels; for the details, please see the document of plot_alpha function in trans_alpha class
    # out$diff_abund  <- t1$plot_diff_abund(plot_type = "ggboxplot", select_taxa = t1$plot_diff_bar_taxa , y_start = 0.05, y_increase = 0.1, color_values = plot_pal)
    
    # out$diff_abund  <- out$diff_abund + scale_y_sqrt() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank()) + ylab("Proportion (sqrt)") 
    # out$combined_plot <- out$diff_bar%>% aplot::insert_right(out$diff_abund, width = 1)
    # out$combined_plot 
    
  }
  if(!require("glmmTMB")) install.packages("glmmTMB")
  
  if ("glmm_beta"  %in% method) #works with taxa_level = all
  {
    t1 <- trans_diff$new(dataset = data, taxa_level = "Species", method = "glmm_beta", p_adjust_method = p_adjust_method,
                         formula = linda_formula , filter_thres = filter_thres) # Time + (1|Subject)" fix_formula
    # View(t1$res_diff)
    
    out$res_diff <- t1$res_diff %<>% 
      group_by(Factors) %>% 
      rstatix::adjust_pvalue(p.col = "P.unadj", method = p_adjust_method, output.col = "P.adj") %>% 
      dplyr::select(-Significance) %>% rstatix::add_significance(p.col = "P.adj", output.col = "Significance")
    
    t1$res_diff %<>% subset(P.adj <= 0.05)
    
    # That's supercool, we can see which taxa are influenced by site or time ... and then how they differ from baseline (TP1)
    # first this approach and then lefse or other test when we see differenfces ...
    
    out$plot_diff_bar <- t1$plot_diff_bar(heatmap_cell = "Estimate", heatmap_sig = "Significance", heatmap_lab_fill = "Coefficient")
    
    out$plot_diff_bar2 <-  t1$plot_diff_bar(color_palette = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")), trans = "log10")
    out$plot_diff_bar3 <-  t1$plot_diff_bar(color_values = c("#053061", "white", "#A50026"), trans = "log10", filter_feature = "", text_y_position = "right", cluster_ggplot = "row")
    
    
    out$plot_diff_bar4 <-  t1$plot_diff_bar(color_values = c("#053061", "white", "#A50026"), trans = "log10", filter_feature = "", text_y_position = "right", cluster_ggplot = "row")
    
    
  }
  
  
  # maaslin2 https://github.com/ChiLiubio/microeco/issues/235
  
  
  
  ####---------------------- return
  
  return(out)
}


#' Phyloseq Classifier Function
#'
#' This function performs classification on a given `phyloseq` object using the `microeco` package. 
#' It includes data preprocessing, feature selection, training, and evaluation of a machine learning model.
#'
#' @param physeq A `phyloseq` object. Default is a subset of `ps_up` where Sample is "Saliva".
#' @param y_response Character. The response variable name. Default is "Time".
#' @param x_predictors Character. The predictor variables. Default is "All".
#' @param prop_train Numeric. Proportion of data to use for training. Default is 3/4.
#' @param method Character. Classification method. Default is "rf" (random forest).
#' @param plot_group Character. Group to plot. Default is "all".
#' @param color_values Named vector. Color palette for plots. Default is `time_pal`.
#' @param ref_train_max_mtry Integer. Maximum number of variables randomly sampled as candidates at each split in the training phase. Default is 5.
#' @param ref_train_ntree Numeric vector. Number of trees in the random forest. Default is c(100, 500, 1000).
#' @param feature_imp_nrep Integer. Number of repetitions for feature importance analysis. Default is 1000.
#' @param boruta_pValue Numeric. p-value threshold for Boruta feature selection. Default is 0.05.
#' @param boruta_maxRuns Integer. Maximum number of Boruta iterations. Default is 300.
#' @param seed Integer. Random seed for reproducibility. Default is 123456.
#' 
#' @return A list containing trained models, confusion matrices, feature importance plots, and prediction results.
#' @import microeco phyloseq file2meco tidyverse doParallel caret randomForest
#' @examples
#' # Example usage
#' result <- phyloseq_classifier(ps_up, y_response = "Group", method = "svmRadial")
#' @export

phyloseq_classifier <- function(physeq = ps_up %>% subset_samples(Sample == "Saliva"),
                                y_response = "Time",
                                x_predictors = "All",
                                prop_train = 3/4,
                                method = "rf", plot_group = "all",
                                color_values = time_pal,
                                ref_train_max_mtry = 5,
                                ref_train_ntree = c(100, 500, 1000),
                                feature_imp_nrep = 1000,
                                boruta_pValue = 0.05,
                                boruta_maxRuns = 300,
                                seed = 123456){

  ####---------------------- Load R package

  require(microeco); require(phyloseq); require(file2meco); require(tidyverse)

  out=NULL

  ####---------------------- Extract data
  # sample_data(physeq)$temp  <- rnorm(nsamples(physeq), mean=22, sd=6)

  physeq %>%
    file2meco::phyloseq2meco(.) -> data

  data$sample_table -> env_data

  ####---------------------- feature - trans_classifier {microeco}

  # initialize: use "genotype" as response variable
  # x.predictors parameter is used to select the taxa; here we use all the taxa data in d1$taxa_abund
  t1 <- trans_classifier$new(dataset = data, y.response = y_response, x.predictors = x_predictors)

  if (!is.null(boruta_pValue))
  {
    set.seed(seed)
    t1$cal_feature_sel(boruta.maxRuns = boruta_maxRuns, boruta.pValue = boruta_pValue)

  }

  set.seed(seed)
  # generate train and test set
  t1$cal_split(prop.train = prop_train)

  # Before training the model, we run the set_trainControl to invoke the trainControl function of caret package to generate the parameters used for training.
  #Here we use the default parameters in trainControl function.
  set.seed(seed)

  t1$set_trainControl()
  # t1$set_trainControl(
  #   method = "repeatedcv",
  #   classProbs = TRUE,
  #   savePredictions = TRUE)

  # use default parameter method = "rf"
  # require(doParallel)
  # library(caret)
  # library(randomForest)
  n <- parallel::detectCores()/2 # experiment!
  cl <- parallel::makeCluster(n)
  doParallel::registerDoParallel(cl)

  set.seed(seed)
  t1$cal_train(method = ifelse(method == "logistic_regression", "rf", method), max.mtry = ref_train_max_mtry, ntree = ref_train_ntree)

  set.seed(seed)
  t1$cal_predict()

  out$res_train <- t1$res_train

  # plot the confusionMatrix to check out the performance

  out$res_confusion_stats <- t1$res_confusion_stats
  out$res_confusion_fit <- t1$res_confusion_fit

  if(method != "logistic_regression")
  {
  out$confusionMatrix <- t1$plot_confusionMatrix()
  }
  # t1$plot_confusion()

  if(method != "logistic_regression")
  {
    t1$cal_ROC(input = "train")
    out$plotROCtrain  <- t1$plot_ROC(plot_type = "ROC", size = 0.5, alpha = 0.7)
    out$plotPRtrain <-  t1$plot_ROC(plot_type = "PR", size = 0.5, alpha = 0.7)
    out$resROCtrain <- t1$res_ROC


  #Using cal_ROC and plot_ROC can get the ROC (Receiver Operator Characteristic) curve.
  # out$Specificitysensitivity() <- t1$res_ROC$res_roc
  # out$RecallPrecision() <- t1$res_ROC$res_pr

  t1$cal_ROC(input = "pred")
  out$plotROCpred  <- t1$plot_ROC(plot_type = "ROC", size = 0.5, alpha = 0.7)
  out$plotPRpred <-  t1$plot_ROC(plot_type = "PR", size = 0.5, alpha = 0.7)
  out$resROCpred <- t1$res_ROC
  }


  # default all groups
  #t1$plot_ROC(size = 0.5, alpha = 0.7)

  # default method in caret package without significance
  # generate significance with rfPermute package

  set.seed(seed)

  if(method != "svmRadial")
  {
  t1$cal_feature_imp(rf_feature_sig = TRUE, num.rep = feature_imp_nrep)
  out$res_feature_imp <- t1$res_feature_imp


  # out$res_feature_imp <-  t1$res_feature_imp()
  if(method != "logistic_regression")
  {
  # default method in caret package without significance
  out$plot_feature_imp1 <- t1$plot_feature_imp(coord_flip = TRUE, colour = "grey5", fill = "grey10", width = 0.6, add_sig = TRUE)
  out$plot_feature_imp2 <- t1$plot_feature_imp(coord_flip = TRUE, colour = "grey5", fill = "grey10", width = 0.6, add_sig = FALSE)

  # rf_sig_show = "MeanDecreaseGini": switch to MeanDecreaseGini
  out$plot_feature_imp3 <- t1$plot_feature_imp(show_sig_group = TRUE, rf_sig_show = "MeanDecreaseGini", coord_flip = TRUE, width = 0.6, add_sig = TRUE, group_aggre = TRUE)
  out$plot_feature_imp4 <- t1$plot_feature_imp(show_sig_group = FALSE, colour = "grey5", fill = "grey10", rf_sig_show = "MeanDecreaseGini", coord_flip = TRUE, width = 0.6, add_sig = FALSE, group_aggre = FALSE)
  }

  # out$res_feature_imp <-  t1$res_feature_imp()
  if(method == "logistic_regression")
  {
    # default method in caret package without significance
    out$plot_feature_imp1 <- t1$plot_feature_imp(coord_flip = TRUE, colour = "grey5", fill = "grey10", width = 0.6, add_sig = TRUE)
    out$plot_feature_imp2 <- t1$plot_feature_imp(coord_flip = TRUE, colour = "grey5", fill = "grey10", width = 0.6, add_sig = FALSE)

    # rf_sig_show = "MeanDecreaseGini": switch to MeanDecreaseGini
    out$plot_feature_imp3 <- t1$plot_feature_imp(show_sig_group = TRUE, rf_sig_show = "IncNodePurity", coord_flip = TRUE, width = 0.6, add_sig = TRUE, group_aggre = TRUE)
    out$plot_feature_imp4 <- t1$plot_feature_imp(show_sig_group = FALSE, colour = "grey5", fill = "grey10", rf_sig_show = "IncNodePurity", coord_flip = TRUE, width = 0.6, add_sig = FALSE, group_aggre = FALSE)
  }
  }
  #
  # # show_sig_group = TRUE: show different colors in groups with different significance labels
  # t1$plot_feature_imp(show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE)
  # t1$plot_feature_imp(show_sig_group = TRUE, coord_flip = TRUE, width = 0.6, add_sig = TRUE)
  # # rf_sig_show = "MeanDecreaseGini": switch to MeanDecreaseGini
  # t1$plot_feature_imp(show_sig_group = TRUE, rf_sig_show = "MeanDecreaseGini", coord_flip = TRUE, width = 0.6, add_sig = TRUE)
  # # group_aggre = FALSE: donot aggregate features for each group
  # t1$plot_feature_imp(show_sig_group = TRUE, rf_sig_show = "MeanDecreaseGini", coord_flip = TRUE, width = 0.6, add_sig = TRUE, group_aggre = TRUE)
  #
  #
  # require Boruta package
  # t1$cal_feature_sel(boruta.maxRuns = boruta_maxRuns, boruta.pValue = 0.01)
  #
  # t2 <- trans_classifier$new(dataset = data, y.response = y_response, x.predictors =  x_predictors)
  # t2$cal_feature_sel(boruta.maxRuns = boruta_maxRuns, boruta.pValue = 0.01)
  # t2$cal_split(prop.train = prop_train)
  # t2$set_trainControl()
  # t2$cal_train()
  # t2$cal_predict()
  # t2$plot_confusionMatrix()
  # t2$cal_ROC()
  # t2$plot_ROC(size = 0.5, alpha = 0.7)


  ####---------------------- return

  out$t1 <- clone(t1)

  return(out)
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note see: https://chiliubio.github.io/microeco_tutorial/diversity-based-class.html
#' @note
#' @note .
#' @return .
#' @export
#' @examples
#'

phyloseq_alpha <- function(physeq, color_groups = NULL, order_x_mean = FALSE, measures = c("Observed", "Shannon", "InvSimpson"), anova_set = NULL,  p_adjust_method = "fdr", group = "SampleType", method = "wilcox"){
  
  ####---------------------- Load R package
  
  require(microeco); require(phyloseq); require(file2meco); require(tidyverse)
  
  ####---------------------- Extract data
  # sample_data(physeq)$temp  <- rnorm(nsamples(physeq), mean=22, sd=6)
  
  physeq %>%
    file2meco::phyloseq2meco(.) -> data
  
  data$sample_table -> env_data
  
  ####---------------------- alpha - trans_alpha {microeco}
  
  t1 <- trans_alpha$new(dataset = data,  group = group, by_group = NULL)
  
  out$data_stat <- t1$data_stat
  out$data_alpha <- t1$data_alpha
  
  ####---------------------- stat
  
  t1$cal_diff(method = method, p_adjust_method = p_adjust_method) #, anova_set = NULL)
  
  out$res <- t1$res_diff
  
  ####---------------------- plot
  
  # t1$plot_alpha(measure = "Chao1", order_x_mean = TRUE, add_sig_text_size = 6)
  
  t1$res_diff %>%
    base::subset(Significance != "ns") -> t1$res_diff
  
  
  if(nrow( t1$res_diff) > 0 )
  {
    
    
    if(is.null(color_groups)) {
      RColorBrewer::brewer.pal(out$data_stat[,group] %>%  unique() %>%  length(), "Dark2") -> color_groups
    }
    
    tmp <- list()
    
    for(i in  out$res %>%  distinct(Measure) %>%  pull() ){
      tmp[[i]] <- t1$plot_alpha(measure = i,
                                # color_values = color_values,
                                add_sig_label = "Significance",
                                add_sig = TRUE,
                                color_values = color_groups,
                                order_x_mean = order_x_mean,
                                add_sig_text_size = 5, xtext_size = 12) +
        theme(plot.margin = unit(c(0.1, 0, 0, 1), "cm"))
    }
    
    out$res_diff_plot <- tmp
    
  }
  
  ####---------------------- return
  
  return(out)
}




#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note see: https://chiliubio.github.io/microeco_tutorial/model-based-class.html#trans_network-class
#' @note
#' @note .
#' @return .
#' @export
#' @examples
#'
#'require(phyloseq); require(tidyverse)
#'data("GlobalPatterns")
#'GlobalPatterns %>%  subset_samples(SampleType %in% c("Feces", "Skin")) %>%  phyloseq_alpha_microeco(group = "SampleType") -> outa
#'data("dataset");dataset %>% file2meco::meco2phyloseq() %>% phyloseq_alpha_microeco(group = "Group",measures = "Chao1") -> outa
#'readRDS((url("https://github.com/fconstancias/DivComAnalyses/blob/9e6ebc69e0c5e136c910186a7136e51af4cecc13/data-raw/ps_invivo.RDS?raw=true" ))) %>% filter_taxa(function(x) sum(x > 0) > 110, TRUE) -> ps

phyloseq_network <- function(physeq, env_cols = c("PAA_ng_ml", "PAG_ng_ml")){
  
  ####---------------------- Load R package
  
  require(microeco); require(phyloseq); require(file2meco); require(tidyverse)
  
  ####---------------------- Extract data
  # sample_data(physeq)$temp  <- rnorm(nsamples(physeq), mean=22, sd=6)
  
  physeq %>%
    file2meco::phyloseq2meco(.) -> data
  
  data$sample_table -> env_data
  
  ####----------------------
  
  t1 <- trans_network$new(dataset = dataset, cor_method = "sparcc", use_sparcc_method = "NetCoMi", filter_thres = 0.001, nThreads = 4, taxa_level = "Species")
  
  
  ####----------------------
  
  t1$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb")
  
  # out$res <- t1$res_diff
  
  ####----------------------
  
  t1$cal_module(
    method = "cluster_fast_greedy",
    module_name_prefix = "M"
  )
  
  ####----------------------
  
  t1$save_network(filepath = "~/Desktop/network.gexf")
  
  ####----------------------
  
  
  t1$cal_network_attr()
  
  t1$res_network_attr -> out$res_network_attr
  
  ####----------------------
  
  t1$get_node_table(node_roles = TRUE)
  
  t1$res_node_table -> out$res_node_table
  
  ####----------------------
  
  t1$get_edge_table()
  
  t1$res_edge_table -> out$res_edge_table
  
  ####----------------------
  
  t1$get_adjacency_matrix(attr = "weight")
  
  t1$res_adjacency_matrix -> out$res_adjacency_matrix
  
  ####----------------------
  
  t1$plot_network(method = "ggraph", node_color = "module")
  
  ####----------------------
  
  t1$cal_eigen()
  
  t1$res_eigen_expla -> out$res_eigen_expla
  
  ####----------------------
  
  t1$plot_taxa_roles(
    use_type = 1,
    roles_color_background = FALSE,
    roles_color_values = NULL,
    add_label = FALSE,
    add_label_group = "Network hubs",
    add_label_text = "name",
    label_text_size = 4,
    label_text_color = "grey50",
    label_text_italic = FALSE,
    plot_module = FALSE,
    x_lim = c(0, 1),
    use_level = "Family",
    show_value = c("z", "p"),
    show_number = 1:10,
    plot_color = "Phylum",
    plot_shape = "taxa_roles",
    plot_size = "Abundance") -> out$plot_taxa_roles_1
  
  ####----------------------
  
  t1$plot_taxa_roles(
    use_type = 2,
    roles_color_background = FALSE,
    roles_color_values = NULL,
    add_label = FALSE,
    add_label_group = "Network hubs",
    add_label_text = "name",
    label_text_size = 4,
    label_text_color = "grey50",
    label_text_italic = FALSE,
    plot_module = FALSE,
    x_lim = c(0, 1),
    use_level = "Class",
    show_value = c("z", "p"),
    show_number = 1:10,
    plot_color = "Class",
    plot_shape = "taxa_roles",
    plot_size = "Abundance") -> out$plot_taxa_roles_2
  
  ####----------------------
  
  # create trans_env object
  t2 <- trans_env$new(dataset = dataset, add_data = dataset$sample_table, env_cols =  env_cols)
  # calculate correlations
  t2$cal_cor(add_abund_table = t1$res_eigen, p_adjust_method = "fdr", by_group = "group", cor_method = "spearman")
  # plot the correlation heatmap
  t2$plot_cor() -> out$plot_cor
  
  t2$res_cor -> out$res_cor
  
  ####---------------------- return
  
  return(out)
}

