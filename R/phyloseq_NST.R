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

phyloseq_NST <- function(physeq,
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

