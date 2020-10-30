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
                                           fittype = "parametric",
                                           contrast = c("diet", "high ALA", "low ALA"),
                                           level_facet = "Class",
                                           Group2 = NULL,
                                           gm_mean = FALSE,
                                           Group_group = Group,
                                           boxplot_colors = NULL)
{
  require(tidyverse);require(DESeq2)
  ps %>%
    speedyseq::tax_glom(taxrank = taxrank) -> ps_temp
  
  taxa_names(ps_temp) <-  tax_table(ps_temp)[,taxrank]
  
  ps_temp %>%
    filter_taxa(function(x) sum(x > sumfilter) > prevfilter*nsamples(ps_temp), TRUE) %>% # keep only ASV with more then 10 reads in more tha 50% of the samples. These are not golden numbers but DESEq can be biaised by ASV abundant in just afew samples. This filtering step is to avoid those biases.
    phyloseq_to_deseq2(as.formula(paste0("~ ", paste(Group, Group2, sep = " ", collapse = " + ")))) -> cds # convert 
  
  
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
  
  
  phyloseq::prune_taxa(da_otus,
                       ps_temp ) %>%
    microbiome::transform("Z") %>%
    plot_heatmap(taxa.label = taxrank,
                 taxa.order = resuls_complete %>% arrange(sign) %>% pull(ASV)      ## ordered according to fold-change
    ) +
    facet_grid(as.formula(paste0(level_facet," ~ ",Group)), scales = "free", space = "free") +
    # scale_fill_gradientn(colours = c("cyan", "black", "red"),
    #                        values = scales::rescale(c(-10, -5, -2, -1, -0.5, -0.05, 0, 0.05, 0.5, 1, 2, 5, 10))) + theme_classic() +
    scale_fill_gradient2(name = "Z-score", low = "#d73027" , mid = "#ffffbf", high = "#1a9850",
                         na.value = "transparent", #trans = scales::log_trans(2),
                         midpoint = 0) + 
    theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) -> heatmap
  
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
  
  out <- list("boxplots"=boxplots,
              "volcano_plot"=volcano_plot,
              "heatmap" =heatmap,
              "results"=resuls_complete)
  
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
                                  padjust.method="fdr", adjustment=3, num.taxa=20, select.variables = cor_variables) -> env.taxa.cor
  
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
  return(list("plot" = p ,
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



# phyloseq_correlate_taxa_old <- function(ps_tmp,
#                                         grouping_column,
#                                         adjustment= 3,
#                                         cor_variables,
#                                         method)
# {
#   ps_tmp %>%
#     tax_table() %>%
#     data.frame() -> tmp
#   
#   ps_tmp %>%
#     transform_sample_counts(function(x) x/sum(x) * 100) %>%
#     microbiome::transform("log10") %>%
#     phyloseq_taxa_env_correlation(grouping_column= grouping_column, method= method, pvalue.threshold=0.05,
#                                   padjust.method="fdr", adjustment=3, num.taxa=20, select.variables = cor_variables) -> env.taxa.cor
#   
#   # plot
#   p <- phyloseq_plot_taxa_env_correlation(env.taxa.cor)
#   
#   p$data %>%
#     dplyr::left_join(tmp %>% rownames_to_column("ASV"),
#                      by = c("Taxa" = "ASV")) %>%
#     dplyr::select(-Taxa) %>%
#     dplyr::rename(Taxa = Strain) %>%
#     phyloseq_plot_taxa_env_correlation() +
#     scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
#                          midpoint = 0, limit = c(-1,1), space = "Lab") -> plot
#   
#   return(list("plot" = plot ,
#               "table" = env.taxa.cor))
# }
# 

