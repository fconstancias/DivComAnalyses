#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note Do not transform it to relative abundance table
#' @note group_var: The name of the group indicator. group_var is required for detecting structural zeros and outliers
#' @note lib_cut: Numeric. Samples with library size less than lib_cut are not included in the analysis.
#' @note zero_cut: Numerical fraction between 0 and 1. Taxa with proportion of zeroes greater than zero_cut are not included in the analysis.
#' @note out_cut Numerical fraction between 0 and 1. For each taxon, observations with proportion of mixture distribution less than out_cut will be detected as outlier zeros; while observations with proportion of mixture distribution greater than 1 - out_cut will be detected as outlier values.
#' @note main_var: The name of the main variable of interest
#' @note neg_lb: Logical. TRUE indicates a taxon would be classified as a structural zero in the corresponding experimental group using its asymptotic lower bound
#' @note https://github.com/FrederickHuangLin/ANCOM
#' @return
#' @export
#' @examples
#' @examples require(phyloseq);require(tidyverse)
#' @examples data("GlobalPatterns")
#' @examples source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_taxa_tests.R")
#' @examples GlobalPatterns %>%
#' @examples  subset_samples(SampleType %in% c("Skin", "Feces")) %>%
#' @examples  filter_taxa(function(x) sum(x > 2000) > 1, TRUE) %>%
#' @examples  phyloseq_ANCOM2(phyloseq = .,
#' @examples                  main_var = "SampleType",
#' @examples                alpha = 0.5)

phyloseq_ANCOM2 <- function(phyloseq,
                            sample_var = "tmp_id",
                            group_var = NULL,
                            out_cut = 0.05,
                            zero_cut = 0.90,
                            lib_cut = 1000,
                            main_var = NULL,
                            p_adj_method = "fdr",
                            alpha = 0.05,
                            adj_formula = NULL,
                            neg_lb = FALSE,
                            rand_formula = NULL,
                            lme_control = NULL){

  #######------------------

  phyloseq %>% otu_table() %>% data.frame() -> feature_table

  phyloseq %>% sample_data() %>% data.frame() %>%
    rownames_to_column(sample_var) -> meta_data

  #######------------------

  # otu_data = read_tsv("~/Documents/GitHub/ANCOM/data/moving-pics-table.tsv", skip = 1)
  # otu_id = otu_data$`feature-id`
  # otu_data = data.frame(otu_data[, -1], check.names = FALSE)
  # rownames(otu_data) = otu_id
  #
  # meta_data = read_tsv("~/Documents/GitHub/ANCOM/data/moving-pics-sample-metadata.tsv")[-1, ]
  # meta_data = meta_data %>%
  #   rename(Sample.ID = SampleID)
  #
  # feature_table = otu_data; sample_var = "Sample.ID"; group_var = NULL
  #
  # main_var = "Subject"; p_adj_method = "BH"; alpha = 0.05
  # adj_formula = NULL; rand_formula = NULL; lme_control = NULL

  #######------------------


  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var,
                                     out_cut, zero_cut, lib_cut, neg_lb)

  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info

  #######------------------

  res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method,
              alpha, adj_formula, rand_formula, lme_control)

  #######------------------

  return(res)

}


library(nlme)
library(dplyr)
library(ggplot2)
library(compositions)

# OTU table should be a matrix/data.frame with each feature in rows and sample in columns.
# Metadata should be a matrix/data.frame containing the sample identifier.

# Data Pre-Processing
feature_table_pre_process = function(feature_table, meta_data, sample_var, group_var = NULL,
                                     out_cut = 0.05, zero_cut = 0.90, lib_cut, neg_lb){
  feature_table = data.frame(feature_table, check.names = FALSE)
  meta_data = data.frame(meta_data, check.names = FALSE)
  # Drop unused levels
  meta_data[] = lapply(meta_data, function(x) if(is.factor(x)) factor(x) else x)
  # Match sample IDs between metadata and feature table
  sample_ID = intersect(meta_data[, sample_var], colnames(feature_table))
  feature_table = feature_table[, sample_ID]
  meta_data = meta_data[match(sample_ID, meta_data[, sample_var]), ]

  # 1. Identify outliers within each taxon
  if (!is.null(group_var)) {
    group = meta_data[, group_var]
    z = feature_table + 1 # Add pseudo-count (1)
    f = log(z)
    f[f == 0] = NA
    f = colMeans(f, na.rm = T)
    f_fit = lm(f ~ group)
    e = rep(0, length(f))
    e[!is.na(group)] = residuals(f_fit)
    y = t(t(z) - e)

    outlier_check = function(x){
      # Fitting the mixture model using the algorithm of Peddada, S. Das, and JT Gene Hwang (2002)
      mu1 = quantile(x, 0.25, na.rm = T)
      mu2 = quantile(x, 0.75, na.rm = T)
      sigma1 = quantile(x, 0.75, na.rm = T) - quantile(x, 0.25, na.rm = T)
      sigma2 = sigma1
      pi = 0.75
      n = length(x)
      epsilon = 100
      tol = 1e-5
      score = pi*dnorm(x, mean = mu1, sd = sigma1)/((1 - pi)*dnorm(x, mean = mu2, sd = sigma2))
      while (epsilon > tol) {
        grp1_ind = (score >= 1)
        mu1_new = mean(x[grp1_ind]); mu2_new = mean(x[!grp1_ind])
        sigma1_new = sd(x[grp1_ind]); if(is.na(sigma1_new)) sigma1_new = 0
        sigma2_new = sd(x[!grp1_ind]); if(is.na(sigma2_new)) sigma2_new = 0
        pi_new = sum(grp1_ind)/n

        para = c(mu1_new, mu2_new, sigma1_new, sigma2_new, pi_new)
        if(any(is.na(para))) break

        score = pi_new * dnorm(x, mean = mu1_new, sd = sigma1_new)/
          ((1-pi_new) * dnorm(x, mean = mu2_new, sd = sigma2_new))

        epsilon = sqrt((mu1 - mu1_new)^2 + (mu2 - mu2_new)^2 +
                         (sigma1 - sigma1_new)^2 + (sigma2 - sigma2_new)^2 + (pi - pi_new)^2)
        mu1 = mu1_new; mu2 = mu2_new; sigma1 = sigma1_new; sigma2 = sigma2_new; pi = pi_new
      }

      if(mu1 + 1.96 * sigma1 < mu2 - 1.96 * sigma2){
        if(pi < out_cut){
          out_ind = grp1_ind
        }else if(pi > 1 - out_cut){
          out_ind = (!grp1_ind)
        }else{
          out_ind = rep(FALSE, n)
        }
      }else{
        out_ind = rep(FALSE, n)
      }
      return(out_ind)
    }
    out_ind = matrix(FALSE, nrow = nrow(feature_table), ncol = ncol(feature_table))
    out_ind[, !is.na(group)] = t(apply(y, 1, function(i)
      unlist(tapply(i, group, function(j) outlier_check(j)))))

    feature_table[out_ind] = NA
  }

  # 2. Discard taxa with zeros  >=  zero_cut
  zero_prop = apply(feature_table, 1, function(x) sum(x == 0, na.rm = T)/length(x[!is.na(x)]))
  taxa_del = which(zero_prop >= zero_cut)
  if(length(taxa_del) > 0){
    feature_table = feature_table[- taxa_del, ]
  }

  # 3. Discard samples with library size < lib_cut
  lib_size = colSums(feature_table, na.rm = T)
  if(any(lib_size < lib_cut)){
    subj_del = which(lib_size < lib_cut)
    feature_table = feature_table[, - subj_del]
    meta_data = meta_data[- subj_del, ]
  }

  # 4. Identify taxa with structure zeros
  if (!is.null(group_var)) {
    group = factor(meta_data[, group_var])
    present_table = as.matrix(feature_table)
    present_table[is.na(present_table)] = 0
    present_table[present_table != 0] = 1

    p_hat = t(apply(present_table, 1, function(x)
      unlist(tapply(x, group, function(y) mean(y, na.rm = T)))))
    samp_size = t(apply(feature_table, 1, function(x)
      unlist(tapply(x, group, function(y) length(y[!is.na(y)])))))
    p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)

    struc_zero = (p_hat == 0) * 1
    # Whether we need to classify a taxon into structural zero by its negative lower bound?
    if(neg_lb) struc_zero[p_hat_lo <= 0] = 1

    # Entries considered to be structural zeros are set to be 0s
    struc_ind = struc_zero[, group]
    feature_table = feature_table * (1 - struc_ind)

    colnames(struc_zero) = paste0("structural_zero (", colnames(struc_zero), ")")
  }else{
    struc_zero = NULL
  }

  # 5. Return results
  res = list(feature_table = feature_table, meta_data = meta_data, structure_zeros = struc_zero)
  return(res)
}

# ANCOM main function
ANCOM = function(feature_table, meta_data, struc_zero = NULL, main_var, p_adj_method = "BH",
                 alpha = 0.05, adj_formula = NULL, rand_formula = NULL, lme_control = NULL){
  # OTU table transformation:
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  if (!is.null(struc_zero)) {
    num_struc_zero = apply(struc_zero, 1, sum)
    comp_table = feature_table[num_struc_zero == 0, ]
  }else{
    comp_table = feature_table
  }
  comp_table = log(as.matrix(comp_table) + 1)
  n_taxa = dim(comp_table)[1]
  taxa_id = rownames(comp_table)
  n_samp = dim(comp_table)[2]

  # Determine the type of statistical test and its formula.
  if (is.null(rand_formula) & is.null(adj_formula)) {
    # Basic model
    # Whether the main variable of interest has two levels or more?
    if (length(unique(meta_data%>%pull(main_var))) == 2) {
      # Two levels: Wilcoxon rank-sum test
      tfun = stats::wilcox.test
    } else{
      # More than two levels: Kruskal-Wallis test
      tfun = stats::kruskal.test
    }
    # Formula
    tformula = formula(paste("x ~", main_var, sep = " "))
  }else if (is.null(rand_formula) & !is.null(adj_formula)) {
    # Model: ANOVA
    tfun = stats::aov
    # Formula
    tformula = formula(paste("x ~", main_var, "+", adj_formula, sep = " "))
  }else if (!is.null(rand_formula)) {
    # Model: Mixed-effects model
    tfun = nlme::lme
    # Formula
    if (is.null(adj_formula)) {
      # Random intercept model
      tformula = formula(paste("x ~", main_var))
    }else {
      # Random coefficients/slope model
      tformula = formula(paste("x ~", main_var, "+", adj_formula))
    }
  }

  # Calculate the p-value for each pairwise comparison of taxa.
  p_data = matrix(NA, nrow = n_taxa, ncol = n_taxa)
  colnames(p_data) = taxa_id
  rownames(p_data) = taxa_id
  pb = txtProgressBar(0, n_taxa - 1, style = 3)
  for (i in 1:(n_taxa - 1)) {
    setTxtProgressBar(pb, i)
    # Loop through each taxon.
    # For each taxon i, additive log ratio (alr) transform the OTU table using taxon i as the reference.
    # e.g. the first alr matrix will be the log abundance data (comp_table) recursively subtracted
    # by the log abundance of 1st taxon (1st column) column-wisely, and remove the first i columns since:
    # the first (i - 1) columns were calculated by previous iterations, and
    # the i^th column contains all zeros.
    alr_data = apply(comp_table, 1, function(x) x - comp_table[i, ])
    # apply(...) allows crossing the data in a number of ways and avoid explicit use of loop constructs.
    # Here, we basically want to iteratively subtract each column of the comp_table by its i^th column.
    alr_data = alr_data[, - (1:i), drop = FALSE]
    n_lr = dim(alr_data)[2] # number of log-ratios (lr)
    alr_data = cbind(alr_data, meta_data) # merge with the metadata

    # P-values
    if (is.null(rand_formula) & is.null(adj_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        test_data = data.frame(x, alr_data, check.names = FALSE)
        suppressWarnings(p <- tfun(tformula, data = test_data)$p.value)
        return(p)
      }
      )
    }else if (is.null(rand_formula) & !is.null(adj_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        fit = tfun(tformula,
                   data = data.frame(x, alr_data, check.names = FALSE),
                   na.action = na.omit)
        p = summary(fit)[[1]][main_var, "Pr(>F)"]
        return(p)
      }
      )
    }else if (!is.null(rand_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        fit = try(tfun(fixed = tformula,
                       data = data.frame(x, alr_data, check.names = FALSE),
                       random = formula(rand_formula),
                       na.action = na.omit,
                       control = lme_control),
                  silent = TRUE)

        if (inherits(fit, "try-error")) {
          p = NA
        } else {
          p = anova(fit)[main_var, "p-value"]
        }
        return(p)
      }
      )
    }
  }
  close(pb)
  # Complete the p-value matrix.
  # What we got from above iterations is a lower triangle matrix of p-values.
  p_data[upper.tri(p_data)] = t(p_data)[upper.tri(p_data)]
  diag(p_data) = 1 # let p-values on diagonal equal to 1
  p_data[is.na(p_data)] = 1 # let p-values of NA equal to 1

  # Multiple comparisons correction.
  q_data = apply(p_data, 2, function(x) p.adjust(x, method = p_adj_method))

  # Calculate the W statistic of ANCOM.
  # For each taxon, count the number of q-values < alpha.
  W = apply(q_data, 2, function(x) sum(x < alpha))

  # Organize outputs
  out_comp = data.frame(taxa_id, W, row.names = NULL, check.names = FALSE)
  # Declare a taxon to be differentially abundant based on the quantile of W statistic.
  # We perform (n_taxa - 1) hypothesis testings on each taxon, so the maximum number of rejections is (n_taxa - 1).
  out_comp = out_comp %>%
    mutate(detected_0.9 = ifelse(W > 0.9 * (n_taxa -1), TRUE, FALSE),
           detected_0.8 = ifelse(W > 0.8 * (n_taxa -1), TRUE, FALSE),
           detected_0.7 = ifelse(W > 0.7 * (n_taxa -1), TRUE, FALSE),
           detected_0.6 = ifelse(W > 0.6 * (n_taxa -1), TRUE, FALSE))

  # Taxa with structural zeros are automatically declared to be differentially abundant
  if (!is.null(struc_zero)){
    out = data.frame(taxa_id = rownames(struc_zero), W = Inf, detected_0.9 = TRUE,
                     detected_0.8 = TRUE, detected_0.7 = TRUE, detected_0.6 = TRUE,
                     row.names = NULL, check.names = FALSE)
    out[match(taxa_id, out$taxa_id), ] = out_comp
  }else{
    out = out_comp
  }

  # Draw volcano plot
  # Calculate clr
  clr_table = apply(feature_table, 2, clr)
  # Calculate clr mean difference
  eff_size = apply(clr_table, 1, function(y)
    lm(y ~ x, data = data.frame(y = y,
                                x = meta_data %>% pull(main_var),
                                check.names = FALSE))$coef[-1])

  if (is.matrix(eff_size)){
    # Data frame for the figure
    dat_fig = data.frame(taxa_id = out$taxa_id, t(eff_size), y = out$W, check.names = FALSE) %>%
      mutate(zero_ind = factor(ifelse(is.infinite(y), "Yes", "No"), levels = c("Yes", "No"))) %>%
      gather(key = group, value = x, rownames(eff_size))
    # Replcace "x" to the name of covariate
    dat_fig$group = sapply(dat_fig$group, function(x) gsub("x", paste0(main_var, " = "), x))
    # Replace Inf by (n_taxa - 1) for structural zeros
    dat_fig$y = replace(dat_fig$y, is.infinite(dat_fig$y), n_taxa - 1)

    fig = ggplot(data = dat_fig) + aes(x = x, y = y) +
      geom_point(aes(color = zero_ind)) +
      facet_wrap(~ group) +
      labs(x = "CLR mean difference", y = "W statistic") +
      scale_color_discrete(name = "Structural zero", drop = FALSE) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "top",
            strip.background = element_rect(fill = "white"))
    fig
  } else{
    # Data frame for the figure
    dat_fig = data.frame(taxa_id = out$taxa_id, x = eff_size, y = out$W) %>%
      mutate(zero_ind = factor(ifelse(is.infinite(y), "Yes", "No"), levels = c("Yes", "No")))
    # Replace Inf by (n_taxa - 1) for structural zeros
    dat_fig$y = replace(dat_fig$y, is.infinite(dat_fig$y), n_taxa - 1)

    fig = ggplot(data = dat_fig) + aes(x = x, y = y) +
      geom_point(aes(color = zero_ind)) +
      labs(x = "CLR mean difference", y = "W statistic") +
      scale_color_discrete(name = "Structural zero", drop = FALSE) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "top")
    fig
  }

  res = list(p_data = p_data, q_data = q_data, out = out, fig = fig)
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
#'

phyloseq_run_DESeq2_pair_plots_formula <- function(ps,
                                                   taxrank = "Strain",
                                                   sumfilter = 10,
                                                   prevfilter = 0.33,
                                                   Group = "diet",
                                                   formula =  paste("age", "deit", sep = " ", collapse = " + "),
                                                   fittype = "parametric",
                                                   taxnames_rm = c("unknown, Incertae Sedis"),
                                                   contrast = c("diet", "high ALA", "low ALA"),
                                                   level_facet = "Class",
                                                   Group2 = NULL,
                                                   gm_mean = FALSE,
                                                   Group_group = Group,
                                                   boxplot_colors = NULL)
{
  require(tidyverse);require(DESeq2)

  ############ ----------------------------
  ps %>%
    filter_taxa(function(x) sum(x > 0) > 0, TRUE)  %>%
    speedyseq::tax_glom(taxrank = taxrank) -> ps_temp

  print(paste0("from ",
               ntaxa(ps)," features, ",
               ntaxa(ps_temp)," were kept after removing features with 0 abundances and  taxa agglomeration"))


  ############ ----------------------------

  if (taxrank != "Strain"){
    prune_taxa(data.frame(tax_table(ps_temp)[,taxrank])  %>%
                 dplyr::filter(!get(taxrank) %in% taxnames_rm) %>% rownames(),ps_temp) -> ps_temp

    taxa_names(ps_temp) <-  tax_table(ps_temp)[,taxrank]
  }


  ############ ----------------------------

  ps_temp %>%
    filter_taxa(function(x){sum(x > sumfilter) >  prevfilter*nsamples(ps_temp)}, prune = TRUE) -> ps_filtered

  ps_filtered  %>%
    phyloseq_to_deseq2(as.formula(paste0("~ ", formula ))) -> cds # convert

  print(paste0("from ", ntaxa(ps_temp)," features, ", ntaxa(ps_filtered)," were kept after taxa agglomeration, sum filter and prevalence filtering"))

  ############ ----------------------------

  if(gm_mean)
  {
    gm_mean = function(x, na.rm = TRUE) {
      exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
    }
    geoMeans = apply(counts(cds), 1, gm_mean)
    cds = estimateSizeFactors(cds, geoMeans = geoMeans)
  }


  ############ ----------------------------run DESeq function

  cds %>%
    DESeq(fitType = fittype) -> dds

  dds %>%
    results(contrast = contrast,
            tidy = TRUE) %>%
    dplyr::rename('ASV' = row) %>%
    left_join(tax_table(ps_temp) %>% data.frame() %>% rownames_to_column('ASV'),
              by = "ASV",
              suffix = c("", ".y")) -> results

  results %>%
    filter(padj < 0.05) %>%
    pull(ASV) -> da_otus

  ############ ----------------------------

  results %>%
    mutate(abs_log2FoldChange = abs(log2FoldChange)) %>% # create new column abs_log2FoldChange absolute values of log2FoldChange column.
    mutate(SIGN  = ifelse(padj <=0.05 & abs_log2FoldChange > 0 , "SIGN", "NS")) %>% # create new column SIGN for each ASV SIGN is added if padj is <=0.05  and abs_log2FoldChange >1.
    mutate(sign = ifelse(log2FoldChange <0, "neg", "pos")) %>% # create new column sign column to specify if  log2FoldChange is positive or negative.
    drop_na(SIGN, # remove NA values for SIGN, log2FoldChange and padj columns.
            log2FoldChange,
            padj) -> resuls_complete

  ############ ----------------------------

  if(length(da_otus)>0)
  {
    phyloseq::prune_taxa(da_otus,
                         ps_temp ) %>%
      microbiome::transform("Z") %>%
      plot_heatmap(taxa.label = taxrank,
                   taxa.order = resuls_complete %>% arrange(sign) %>% pull(ASV),
                   method = "MDS", distance = "euclidean") +
      facet_grid(as.formula(paste0(level_facet," ~ ",Group)), scales = "free", space = "free") +
      # scale_fill_gradientn(colours = c("cyan", "black", "red"),
      #                        values = scales::rescale(c(-10, -5, -2, -1, -0.5, -0.05, 0, 0.05, 0.5, 1, 2, 5, 10))) + theme_classic() +
      scale_fill_gradient2(name = "Z-score", low = "#d73027" , mid = "#ffffbf", high = "#1a9850",
                           na.value = "transparent", #trans = scales::log_trans(2),
                           midpoint = 0) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) -> heatmap

    ############ ----------------------------

    phyloseq::prune_taxa(da_otus,
                         ps_temp ) %>%
      transform_sample_counts(function(x) x/sum(x) * 100) %>%
      plot_heatmap(taxa.label = taxrank,
                   taxa.order = resuls_complete %>% arrange(sign) %>% pull(ASV),
                   method = "MDS", distance = "euclidean") +
      facet_grid(as.formula(paste0(level_facet," ~ ",Group)), scales = "free", space = "free") +
      # scale_fill_gradientn(colours = c("cyan", "black", "red"),
      #                        values = scales::rescale(c(-10, -5, -2, -1, -0.5, -0.05, 0, 0.05, 0.5, 1, 2, 5, 10))) + theme_classic() +
      scale_fill_gradient(name = "Proportion - %", low = "#d73027" , mid = "#ffffbf", high = "#1a9850",
                          na.value = "transparent") +
      theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) -> heatmap_prop

    ############ ----------------------------

    resuls_complete %>%
      ggplot(aes(x = log2FoldChange, y = -log10(padj))) + # tell ggplot that we are going to plot -log10(padj) as a function of log2FoldChange
      geom_point(aes(shape = SIGN, # points are foing to be ploted, shape is coded by SIGN column (SIGN or NS)
                     size = SIGN, # size is coded by SIGN column (SIGN or NS)
                     fill = get(level_facet), # filling colour and line color is coded by Class info of the ASV
                     colour = get(level_facet),
                     alpha = baseMean)) + # alpha = transparency reflects baseMean i.e., mean value of ASV among samples.
      scale_shape_manual(values = c(4, 21)) + # We force the shape of the points to be 4 and 21.see : <http://www.sthda.com/sthda/RDoc/images/points-symbols.png>
      scale_alpha_continuous(name = "baseMean",
                             limits = c(0,1000),
                             trans = "sqrt",
                             range = c(0.6, 0.8)) + #  transparency values from 06 to 0.8
      scale_colour_viridis_d(alpha = 0.7,
                             begin = 0,
                             end = 1,
                             direction = 1) +
      scale_fill_viridis_d() +
      scale_size_manual(values=c(0.2, 2)) + # We force the size of the points: 0.2 for NS and 2 for SIGN
      # geom_text_repel( # We use ggrepel to display Strain column for significant ASVs
      #   data = resuls_complete %>%
      #     drop_na(SIGN, log2FoldChange ,padj) %>%
      #     subset(padj <= 0.001 & abs_log2FoldChange > 4),
      #   aes(label = Genus),
      #   size = 2,
      #   force = 4,
      # ) +
      geom_hline( # adding horizontal line:
        yintercept = -log10(0.05),
        col = "red",
        linetype = "dotted",
        size = 0.5
      ) + geom_vline( # adding vertical lines:
        xintercept = c(-2, 2),
        col = "red",
        linetype = "dotted",
        size = 0.5) -> volcano_plot


    prune_taxa(da_otus,
               ps_temp %>% transform_sample_counts(function(x) x/sum(x) * 100)) -> ps_tmp #%>%
    # subset_taxa(Family != "unknown")-> ps_tmp

    taxa_names(ps_tmp) <- tax_table(ps_tmp)[, taxrank]

    lapply(
      as.list(taxa_names(ps_tmp)),
      FUN = phyloseq_boxplot_abundance,
      ps = ps_tmp,
      x= Group, color = Group_group, level = taxrank, line=NULL, violin = FALSE, show.points = TRUE, colors = boxplot_colors) -> boxplots

    names(boxplots) <- taxa_names(ps_tmp)

    out <- list("ps_filtered" = ps_filtered,
                "boxplots"=boxplots,
                "volcano_plot"=volcano_plot,
                "heatmap" =heatmap,
                "heatmap_prop" = heatmap_prop,
                "results"=resuls_complete)

  }else{
    print("No singinifcant features found")

    out <- list("ps_filtered" = ps_filtered,
                "results"= resuls_complete)
  }
  return(out)
  detach("package:DESeq2", unload = TRUE)

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
#'data("GlobalPatterns")
#'GlobalPatterns %>%
#' subset_samples(SampleType %in% c("Feces", "Skin")) %>%
#' phyloseq_get_strains_fast  -> ps_tmp
#'
#'ps_tmp %>% phyloseq_Maaslin2(fixed_effects = "SampleType",random_effects = NULL , min_abundance = 2, min_prevalence = 0.5 ,  normalization = "NONE", transform = "NONE", analysis_method = "ZINB", output_dir = "~/test_masslin2_ZINB/")
#'ps_tmp %>% phyloseq_Maaslin2(taxrank = FALSE, rename_ASV_strain = FALSE, fixed_effects = "SampleType",random_effects = NULL , min_abundance = 2, min_prevalence = 0.5 ,  normalization = "NONE", transform = "NONE", analysis_method = "NEGBIN", output_dir = "~/test_masslin2_negbin/") -> test
#'ps_tmp %>% phyloseq_Maaslin2(fixed_effects = "SampleType",random_effects = NULL , min_abundance = 2, min_prevalence = 0.5 ,  normalization = "NONE", transform = "NONE", analysis_method = "CPLM", output_dir = "~/test_masslin2_CPLM/") -> test_CPLM


phyloseq_Maaslin2 <- function(phyloseq,
                              min_abundance = 0,
                              min_prevalence = 0.1 ,
                              min_variance = 0,
                              random_effects = NULL,
                              fixed_effects = c("treatment"),
                              max_significance = 0.25,
                              normalization = "TSS",
                              transform = "LOG",
                              analysis_method = "LM",
                              correction = "BH",
                              standardize = TRUE,
                              reference = NULL,
                              cores = 4,
                              plot_heatmap = FALSE,
                              plot_scatter = TRUE,
                              heatmap_first_n = 50,
                              output_dir = "~/test_masslin2/",
                              add_ASV_taxonomy = TRUE){

  ##---------------------------------------------
  require(tidyverse); require(Maaslin2); require(phyloseq)

  ##---------------------------------------------

  Maaslin2(phyloseq %>% otu_table() %>%  t(),
           phyloseq %>% sample_data() %>% data.frame(),
           output_dir,
           analysis_method = analysis_method,
           normalization = normalization,
           transform = transform,
           min_abundance = min_abundance,
           min_prevalence = min_prevalence,
           min_variance = min_variance,
           random_effects = random_effects,
           fixed_effects = fixed_effects,
           correction =  correction,
           reference = reference,
           standardize = standardize,
           max_significance = max_significance,
           cores = cores,
           heatmap_first_n = heatmap_first_n,
           plot_heatmap = plot_heatmap,
           plot_scatter = plot_scatter) -> out

  if(add_ASV_taxonomy == TRUE){

    phyloseq %>%
      tax_table() %>%
      as.data.frame()%>%
      rownames_to_column(var = "feature") -> tax_table

    left_join(out$results,
              tax_table, by="feature") -> out$results_ASV_tax

    write_tsv(out$results_ASV_tax, file=paste0(output_dir, "results_tax_info.tsv"))

  }
  ##---------------------------------------------

  # gc()

  ##---------------------------------------------

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
#'

phyloseq_run_DESeq2_pair_plots <- function(ps,
                                           taxrank = "Strain",
                                           sumfilter = 10,
                                           prevfilter = 0.33,
                                           Group = "diet",
                                           taxnames_rm = c("unknown, Incertae Sedis"),
                                           fittype = "parametric",
                                           contrast = c("diet", "high ALA", "low ALA"),
                                           level_facet = "Class",
                                           Group2 = NULL,
                                           gm_mean = FALSE,
                                           Group_group = Group,
                                           boxplot_colors = NULL,
                                           generate_plots = TRUE)
{
  require(tidyverse);require(DESeq2)
  ps %>%
    filter_taxa(function(x) sum(x > 0) > 0, TRUE)  %>%
    speedyseq::tax_glom(taxrank = taxrank) -> ps_temp

  print(paste0("from ", ntaxa(ps)," features, ", ntaxa(ps_temp)," were kept after removing features with 0 abundances and  taxa agglomeration"))



  if (taxrank != "Strain"){
    prune_taxa(data.frame(tax_table(ps_temp)[,taxrank])  %>%
                 dplyr::filter(!get(taxrank) %in% taxnames_rm) %>% rownames(),ps_temp) -> ps_temp

    taxa_names(ps_temp) <-  tax_table(ps_temp)[,taxrank]
  }


  # ps_temp %>%
  #   microbiome::core(detection = sumfilter, prevalence = prevfilter) -> ps_filtered

  ps_temp %>%
    filter_taxa(function(x){sum(x > sumfilter) >  prevfilter*nsamples(ps_temp)}, prune = TRUE) -> ps_filtered


  ps_filtered %>%
    phyloseq_to_deseq2(as.formula(paste0("~ ", paste(Group, Group2, sep = " ", collapse = " + ")))) -> cds # convert


  print(paste0("from ", ntaxa(ps_temp)," features, ", ntaxa(ps_filtered)," were kept after taxa agglomeration, sum filter and prevalence filtering"))

  if(gm_mean)
  {
    gm_mean = function(x, na.rm = TRUE) {
      exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
    }
    geoMeans = apply(counts(cds), 1, gm_mean)
    cds = estimateSizeFactors(cds, geoMeans = geoMeans)
  }


  # run DESeq function
  cds %>%
    DESeq(fitType = fittype) -> dds

  dds %>%
    results(contrast = contrast,
            tidy = TRUE) %>%
    dplyr::rename('ASV' = row) %>%
    left_join(tax_table(ps_temp) %>% data.frame() %>% rownames_to_column('ASV'),
              by = "ASV",
              suffix = c("", ".y")) -> results

  results %>%
    filter(padj < 0.05) %>%
    pull(ASV) -> da_otus

  results %>%
    mutate(abs_log2FoldChange = abs(log2FoldChange)) %>% # create new column abs_log2FoldChange absolute values of log2FoldChange column.
    mutate(SIGN  = ifelse(padj <=0.05 & abs_log2FoldChange > 0 , "SIGN", "NS")) %>% # create new column SIGN for each ASV SIGN is added if padj is <=0.05  and abs_log2FoldChange >1.
    mutate(sign = ifelse(log2FoldChange <0, "neg", "pos")) %>% # create new column sign column to specify if  log2FoldChange is positive or negative.
    drop_na(SIGN, # remove NA values for SIGN, log2FoldChange and padj columns.
            log2FoldChange,
            padj) -> resuls_complete


  if(length(da_otus) > 0 & generate_plots == TRUE)
  {
    phyloseq::prune_taxa(da_otus,
                         ps_temp ) %>%
      microbiome::transform("Z") %>%
      plot_heatmap(taxa.label = taxrank,
                   taxa.order = resuls_complete %>% arrange(sign) %>% pull(ASV),
                   method = "MDS", distance = "euclidean") +
      facet_grid(as.formula(paste0(level_facet," ~ ",Group)), scales = "free", space = "free") +
      # scale_fill_gradientn(colours = c("cyan", "black", "red"),
      #                        values = scales::rescale(c(-10, -5, -2, -1, -0.5, -0.05, 0, 0.05, 0.5, 1, 2, 5, 10))) + theme_classic() +
      scale_fill_gradient2(name = "Z-score", low = "#d73027" , mid = "#ffffbf", high = "#1a9850",
                           na.value = "transparent", #trans = scales::log_trans(2),
                           midpoint = 0) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) -> heatmap


    phyloseq::prune_taxa(da_otus,
                         ps_temp ) %>%
      transform_sample_counts(function(x) x/sum(x) * 100) %>%
      plot_heatmap(taxa.label = taxrank,
                   taxa.order = resuls_complete %>% arrange(sign) %>% pull(ASV),
                   method = "MDS", distance = "euclidean") +
      facet_grid(as.formula(paste0(level_facet," ~ ",Group)), scales = "free", space = "free") +
      # scale_fill_gradientn(colours = c("cyan", "black", "red"),
      #                        values = scales::rescale(c(-10, -5, -2, -1, -0.5, -0.05, 0, 0.05, 0.5, 1, 2, 5, 10))) + theme_classic() +
      scale_fill_gradient(name = "Proportion - %", low = "#d73027" , high = "#1a9850",
                          na.value = "transparent") +
      theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) -> heatmap_prop

    resuls_complete %>%
      ggplot(aes(x = log2FoldChange, y = -log10(padj))) + # tell ggplot that we are going to plot -log10(padj) as a function of log2FoldChange
      geom_point(aes(shape = SIGN, # points are foing to be ploted, shape is coded by SIGN column (SIGN or NS)
                     size = SIGN, # size is coded by SIGN column (SIGN or NS)
                     fill = get(level_facet), # filling colour and line color is coded by Class info of the ASV
                     colour = get(level_facet),
                     alpha = baseMean)) + # alpha = transparency reflects baseMean i.e., mean value of ASV among samples.
      scale_shape_manual(values = c(4, 21)) + # We force the shape of the points to be 4 and 21.see : <http://www.sthda.com/sthda/RDoc/images/points-symbols.png>
      scale_alpha_continuous(name = "baseMean",
                             limits = c(0,1000),
                             trans = "sqrt",
                             range = c(0.6, 0.8)) + #  transparency values from 06 to 0.8
      scale_colour_viridis_d(alpha = 0.7,
                             begin = 0,
                             end = 1,
                             direction = 1) +
      scale_fill_viridis_d() +
      scale_size_manual(values=c(0.2, 2)) + # We force the size of the points: 0.2 for NS and 2 for SIGN
      # geom_text_repel( # We use ggrepel to display Strain column for significant ASVs
      #   data = resuls_complete %>%
      #     drop_na(SIGN, log2FoldChange ,padj) %>%
      #     subset(padj <= 0.001 & abs_log2FoldChange > 4),
      #   aes(label = Genus),
      #   size = 2,
      #   force = 4,
      # ) +
      geom_hline( # adding horizontal line:
        yintercept = -log10(0.05),
        col = "red",
        linetype = "dotted",
        size = 0.5
      ) + geom_vline( # adding vertical lines:
        xintercept = c(-2, 2),
        col = "red",
        linetype = "dotted",
        size = 0.5) -> volcano_plot


    prune_taxa(da_otus,
               ps_temp %>% transform_sample_counts(function(x) x/sum(x) * 100)) -> ps_tmp #%>%
    # subset_taxa(Family != "unknown")-> ps_tmp

    taxa_names(ps_tmp) <- tax_table(ps_tmp)[, taxrank]

    lapply(
      as.list(taxa_names(ps_tmp)),
      FUN = phyloseq_boxplot_abundance,
      ps = ps_tmp,
      x= Group, color = Group_group, level = taxrank, line=NULL, violin = FALSE, show.points = TRUE, colors = boxplot_colors) -> boxplots

    names(boxplots) <- taxa_names(ps_tmp)

    out <- list("ps_filtered" = ps_filtered,
                "boxplots"=boxplots,
                "volcano_plot"=volcano_plot,
                "heatmap" =heatmap,
                "heatmap_prop" = heatmap_prop,
                "results"=resuls_complete)

  }else{
    print("No singinifcant features found or you did not want any plots?")

    out <- list("ps_filtered" = ps_filtered,
                "results"=resuls_complete)
  }
  return(out)
  detach("package:DESeq2", unload = TRUE)

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


phyloseq_heatmap_boxplots <- function(physeq_mIMT1,
                                      da_otu,
                                      groups = c("mIMT_Dysbiotic", "mIMT_Eubiotic"),
                                      comp = "DysvsEub",
                                      var = "treatment",
                                      taxrank = "Strain",
                                      level_facet = "Class",
                                      trans = "Z",
                                      Group_group = var,
                                      boxplot_colors = NULL){

  ####-------- Extract sample belonging to groups of var

  prune_samples(get_variable(physeq_mIMT1, var) %in% groups,
                physeq_mIMT1) -> ps

  ####-------- Continue if there as significant features

  if(length(da_otu)>0)
  {

    ####-------- transform phyloseq object before selecting features:
    ps %>%
      microbiome::transform(transform = trans) -> ps

    ####-------- generate heatmap of those

    phyloseq::prune_taxa(da_otu,
                         ps) %>%
      plot_heatmap(taxa.label = taxrank, method = "MDS", distance = "euclidean") +
      facet_grid(as.formula(paste0(level_facet," ~ ",var)), scales = "free", space = "free") +
      theme_classic() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7)) + xlab(NULL) + ylab(NULL) -> heatmap


    if(trans == "Z"){
      heatmap + scale_fill_gradient2(name = "Z-score", low = "#d73027" , mid = "#ffffbf", high = "#1a9850",
                                     na.value = "transparent", #trans = scales::log_trans(2),
                                     midpoint = 0) -> heatmap
    }
    if(trans == "compositional"){
      heatmap + scale_fill_gradient(name = "Proportion - %", low = "#d73027" , high = "#1a9850",
                                    na.value = "transparent")  -> heatmap
    }

    ####-------- generate boxplots

    prune_taxa(da_otu,
               ps %>% transform_sample_counts(function(x) x/sum(x) * 100)) -> ps_tmp #%>%
    # subset_taxa(Family != "unknown")-> ps_tmp

    taxa_names(ps_tmp) <- tax_table(ps_tmp)[, taxrank]

    lapply(
      as.list(taxa_names(ps_tmp)),
      FUN = phyloseq_boxplot_abundance,
      ps = ps_tmp,
      x= var, color = Group_group, level = taxrank, line=NULL, violin = FALSE, show.points = TRUE, colors = boxplot_colors) -> boxplots

    names(boxplots) <- taxa_names(ps_tmp)

  }else{
    print("No singinifcant features found")
  }

  out <- list("heatmap" = heatmap,
              "boxplots" = boxplots)
}

#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note columns must be padj, log2FoldChange, baseMean, SIGN
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'

plot_volcano <- function(resuls_complete,
                         level_facet = "Class"){
  ####-------- generate volcano plot

  resuls_complete %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) + # tell ggplot that we are going to plot -log10(padj) as a function of log2FoldChange
  geom_point(aes(shape = SIGN, # points are foing to be ploted, shape is coded by SIGN column (SIGN or NS)
                 size = SIGN, # size is coded by SIGN column (SIGN or NS)
                 fill = get(level_facet), # filling colour and line color is coded by Class info of the ASV
                 colour = get(level_facet),
                 alpha = baseMean)) + # alpha = transparency reflects baseMean i.e., mean value of ASV among samples.
  scale_shape_manual(values = c(4, 21)) + # We force the shape of the points to be 4 and 21.see : <http://www.sthda.com/sthda/RDoc/images/points-symbols.png>
  scale_alpha_continuous(name = "baseMean",
                         limits = c(0,1000),
                         trans = "sqrt",
                         range = c(0.6, 0.8)) + #  transparency values from 06 to 0.8
  scale_colour_viridis_d(alpha = 0.7,
                         begin = 0,
                         end = 1,
                         direction = 1) +
  scale_fill_viridis_d() +
  scale_size_manual(values=c(0.2, 2)) + # We force the size of the points: 0.2 for NS and 2 for SIGN
  # geom_text_repel( # We use ggrepel to display Strain column for significant ASVs
  #   data = resuls_complete %>%
  #     drop_na(SIGN, log2FoldChange ,padj) %>%
  #     subset(padj <= 0.001 & abs_log2FoldChange > 4),
  #   aes(label = Genus),
  #   size = 2,
  #   force = 4,
  # ) +
  geom_hline( # adding horizontal line:
    yintercept = -log10(0.05),
    col = "red",
    linetype = "dotted",
    size = 0.5
  ) + geom_vline( # adding vertical lines:
    xintercept = c(-2, 2),
    col = "red",
    linetype = "dotted",
    size = 0.5) -> volcano_plot

return(volcano_plot)
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
#'library(phyloseq);library(tidyverse)
#'data(GlobalPatterns)
#'GlobalPatterns %>% phyloseq_get_strains() -> ps2
#'ps2 %>% phyloseq_boxplot_abundance(color = NULL)
#'
#'phyloseq_boxplot_abundance(ps = ps2)
#'
phyloseq_boxplot_abundance <- function (ps,
                                        x,
                                        y,
                                        level = "Strain",
                                        color = NULL, # NULL
                                        line = NULL,
                                        violin = FALSE,
                                        na.rm = FALSE,
                                        show.points = TRUE,
                                        size = 2,
                                        alpha = 1,
                                        log10 = TRUE,
                                        colors = NULL)

{
  require(tidyverse); require(microbiome)
  colors = colors
  change <- xvar <- yvar <- linevar <- colorvar <- NULL
  pseq <- ps
  taxa_names(pseq) <- tax_table(pseq)[,level]

  otu <- microbiome::abundances(pseq)
  df <- microbiome::meta(pseq)
  df$xvar <- df[[x]]
  if (!is.factor(df[[x]])) {
    df$xvar <- factor(as.character(df$xvar))
  }
  if (y %in% taxa_names(pseq)) {
    df$yvar <- as.vector(unlist(otu[y, ]))
  }else {
    df$yvar <- as.vector(unlist(sample_data(pseq)[, y]))
  }
  if (na.rm) {
    df <- subset(df, !is.na(xvar))
    df <- subset(df, !is.na(yvar))
  }
  if (nrow(df) == 0) {
    warning("No sufficient data for plotting available. \n            Returning an empty plot.")
    return(ggplot())
  }
  df$xvar <- factor(df$xvar)

  y %>%
    stringr::str_replace("unknown", "un.") %>%
    stringr::str_trunc(35,  side ="center") -> ylab

  p <- ggplot(df, aes(x = xvar, y = yvar)) + theme_classic() + ylab(ylab)
  if (show.points) {
    p <- p + geom_jitter(size = size,
                         alpha = alpha,
                         aes_string(colour = color,
                                    fill = color))
  }
  if (!violin) {
    p <- p + geom_boxplot(outlier.size = 0,
                          outlier.shape = "",
                          aes_string(color = color,
                                     fill = color),
                          alpha = 0.8)
  }else {
    p <- p + geom_violin(fill = NA,
                         aes_string(color = color,
                                    fill = color))
  }
  if (!is.null(line)) {
    df$linevar <- factor(df[[line]])
    df2 <- suppressWarnings(df %>% arrange(linevar, xvar) %>%
                              group_by(linevar) %>% summarise(change = diff(yvar)))
    df$change <- df2$change[match(df$linevar, df2$linevar)]
    df$change <- sign(df$change)
    p <- p + geom_line(data = df, aes(group = linevar, color = change),
                       size = 1) + scale_colour_gradient2(low = "blue",
                                                          mid = "black", high = "red", midpoint = 0, na.value = "grey50",
                                                          guide = "none")
  }

  if (log10) {
    p <- p + scale_y_log10() + ylab(paste0(ylab)) + theme(legend.position = "none")
  }
  if (is.null(colors)) {
    p <- p + xlab(x) + theme(legend.position = "none")
    return(p)
  }else {
    p <- p + xlab(NULL) + theme(legend.position = "none") +
      scale_colour_manual(values = colors) +
      scale_fill_manual(values = colors)
    return(p)
  }


  detach("package:microbiome", unload = TRUE)
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
#'library(phyloseq)
#'data(GlobalPatterns)
#'phyloseq_A_B_ratio()
#'


phyloseq_A_B_ratio <- function(ps = GlobalPatterns,
                               level = "Phylum",
                               a_name = "Bacteroidetes",
                               b_name = "Firmicutes",
                               Group = "SampleType",
                               sampleID = FALSE)
{
  require(microbiome)
  require(tidyverse)
  require(ggpubr)

  microbiome::transform(microbiome::aggregate_taxa(ps, level = level), "compositional") -> tmp
  a <- microbiome::abundances(tmp)[a_name, ]
  b <- microbiome::abundances(tmp)[b_name, ]

  name <- paste0(a_name, "_", b_name)

  a/b %>%
    data.frame() %>%
    dplyr::rename(!!name :=  ".") -> tmp

  tmp %>%
    tibble::rownames_to_column('sample') %>%
    dplyr::full_join(
      ps %>%
        sample_data() %>%
        data.frame() %>%
        rownames_to_column('sample'),
      by = c("sample" = "sample")
    ) -> df


  df %>%
    ggplot(aes_string(x=Group,
                      y=name,
                      colour=Group,
                      fill = Group,
                      label = "sample")) +
    # facet_grid(as.formula(paste("~","diet")), drop=T,scale="free",space="free_x") +
    geom_boxplot(outlier.colour = NA,alpha=0.8,
                 position = position_dodge(width=0.7)) +
    # geom_violin(alpha = 0.1) +
    geom_jitter(size=2, alpha=1, position=position_jitterdodge(1)) +
    ylab(paste0(paste0(a_name, "/", b_name) , " ratio"))  + xlab(NULL)  +
    theme(axis.text.x = element_blank()) +
    theme_classic() -> p

  if(sampleID == TRUE)
  {
    p +
      ggrepel::geom_text_repel(position = position_jitterdodge(1),
                               size = 2,
                               # color = 'black',
                               segment.color = 'grey50'# ,    min.segment.length = 0
      ) -> p
  }

  ggpubr::compare_means(as.formula(paste0(name, " ~ ", Group)),
                        # group.by = "variable",
                        data = df,
                        method = "wilcox.test") -> KW_tests

  return(out = list("plot" = p,
                    "df" = df,
                    "KW_tests" = KW_tests))

  detach("package:microbiome", unload = TRUE); detach("package:ggpubr", unload = TRUE)
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
#'enterotype %>%
#' subset_samples(!is.na(Age)) %>%
#' phyloseq_correlate_taxa(log10 = TRUE, tax_glom = "Genus", grouping_column = "Gender", cor_variables = "Age")
#'

phyloseq_correlate_taxa <- function(ps_tmp,
                                    log10 = TRUE,
                                    tax_glom = FALSE,
                                    grouping_column,
                                    adjustment= 3,
                                    num_taxa = 20,
                                    cor_variables,
                                    method = "spearman")
{
  require(tidyverse)

  if(tax_glom!=FALSE)
  {
    ps_tmp %>%
      tax_glom(taxrank = tax_glom) -> ps_tmp
  }

  ps_tmp %>%
    transform_sample_counts(function(x) x/sum(x) * 100) -> tmp2

  if(log10==TRUE)
  {
    tmp2 %>%
      microbiome::transform("log10") -> tmp2
  }
  tmp2 %>%
    phyloseq_taxa_env_correlation(grouping_column= grouping_column, method= method, pvalue.threshold=0.05,
                                  padjust.method="fdr", adjustment=adjustment, num.taxa= num_taxa, select.variables = cor_variables) -> env.taxa.cor

  # plot
  p <- phyloseq_plot_taxa_env_correlation(env.taxa.cor)

  if(tax_glom==FALSE)
  {
    as(tax_table(ps_tmp), "matrix") %>%
      data.frame() -> tmp

    p$data %>%
      dplyr::left_join(tmp %>% rownames_to_column("ASV"),
                       by = c("Taxa" = "ASV")) %>%
      dplyr::select(-Taxa) %>%
      dplyr::rename(Taxa = Strain) %>%
      phyloseq_plot_taxa_env_correlation() +
      scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                           midpoint = 0, limit = c(-1,1), space = "Lab") -> p
  }

  if(tax_glom!=FALSE)
  {
    as(tax_table(ps_tmp), "matrix") %>%
      data.frame() -> tmp

    p$data %>%
      dplyr::left_join(tmp %>% rownames_to_column("ASV"),
                       by = c("Taxa" = "ASV")) %>%
      dplyr::select(-Taxa) %>%
      dplyr::rename(Taxa := !!tax_glom) %>%
      phyloseq_plot_taxa_env_correlation() +
      scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                           midpoint = 0, limit = c(-1,1), space = "Lab") -> p
  }

  return(list("plot" = p,
              "table" = env.taxa.cor))
}


phyloseq_plot_taxa_env_correlation <- function(df)
{
  p <- ggplot2::ggplot(aes(x = Type, y = Taxa, fill = Correlation),
                       data = df)
  p <- p + ggplot2::geom_tile() + scale_fill_gradient2(low = "#2C7BB6",
                                                       mid = "white", high = "#D7191C")
  p <- p + ggplot2::theme(axis.text.x = element_text(angle = 90,
                                                     hjust = 1, vjust = 0.5))
  p <- p + ggplot2::geom_text(aes(label = Significance), color = "black",
                              size = 3) + labs(y = NULL, x = NULL)
  p <- p + ggplot2::facet_grid(. ~ Env, drop = TRUE, scale = "free",
                               space = "free_x")
  p <- p + ggplot2::xlab("Groups")
  p <- p + ggplot2::theme(strip.background = element_rect(fill = "white"))
  return(p)
}


phyloseq_taxa_env_correlation <- function (physeq, grouping_column, method = "pearson", pvalue.threshold = 0.05,
                                           padjust.method = "fdr", adjustment = 3, num.taxa = 50, select.variables = NULL)
{
  require(tidyverse)
  method <- match.arg(method, c("pearson", "kendall", "spearman"),
                      several.ok = F)
  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  abund_table <- as.data.frame(otu_table(physeq))
  meta_table <- data.frame(sample_data(physeq))
  groups <- meta_table[, grouping_column]
  if (!is.null(select.variables)) {
    meta_table <- subset(meta_table, select = select.variables)
  }
  # mt_env <- meta_table[, sapply(meta_table, as.numeric)] %>%
  #   column_to_rownames("tmp") %>%
  #   rownames_to_column("tmp")
  mt_env <- meta_table
  abund_table_filt <- abund_table[rownames(mt_env), ]
  abund_table_filt <- abund_table_filt[, order(colSums(abund_table_filt),
                                               decreasing = TRUE)]
  taxa_list <- colnames(abund_table_filt)[1:num.taxa]
  taxa_list <- taxa_list[!grepl("Unknown", taxa_list)]
  abund_table_filt <- data.frame(abund_table_filt[, colnames(abund_table_filt) %in%
                                                    taxa_list])
  df <- tables.correlate(abund_table_filt, mt_env, groups,
                         method)
  colnames(df) <- c("Taxa", "Env", "Correlation", "Pvalue",
                    "Type")
  df$Pvalue <- as.numeric(as.character(df$Pvalue))
  df$Correlation <- as.numeric(as.character(df$Correlation))
  df$AdjPvalue <- rep(0, dim(df)[1])
  df <- p.adjust.cor(df, adjustment, padjust.method)
  df$Significance <- cut(df$AdjPvalue, breaks = c(-Inf, 0.001,
                                                  0.01, 0.05, Inf), label = c("***", "**", "*", ""))
  df <- df[complete.cases(df), ]
  return(df)
}

tables.correlate<-function(table1, table2, groups=NULL, method){
  df<-NULL
  for(i in colnames(table1)){
    for(j in colnames(table2)){

      if(!is.null(groups)){
        for(k in unique(groups)){
          a<-table1[groups==k,i,drop=F]
          b<-table2[groups==k,j,drop=F]
          tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value,k)

          if(is.null(df)){df<-tmp} else{df<-rbind(df,tmp)}
        }
      }
      else{

        a<-table1[,i,drop=F]
        b<-table2[,j,drop=F]
        tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value)

        if(is.null(df)){df<-tmp} else{df<-rbind(df,tmp)}

      }

    }
  }

  df<-data.frame(row.names=NULL,df)
  return(df)
}

# df is a data frame
p.adjust.cor <- function(df,adjustment=1,padjust.method="BH"){
  if(adjustment==1){
    df$AdjPvalue<-df$Pvalue
  } else if (adjustment==2){
    for(i in unique(df$Env)){
      for(j in unique(df$Type)){
        sel<-df$Env==i & df$Type==j
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method=padjust.method)
      }
    }
  } else if (adjustment==3){
    for(i in unique(df$Taxa)){
      for(j in unique(df$Type)){
        sel<-df$Taxa==i & df$Type==j
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method=padjust.method)
      }
    }
  } else if (adjustment==4){
    for(i in unique(df$Taxa)){
      sel<-df$Taxa==i
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method=padjust.method)
    }
  } else if (adjustment==5){
    for(i in unique(df$Env)){
      sel<-df$Env==i
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method=padjust.method)
    }
  }
  return(df)
}


