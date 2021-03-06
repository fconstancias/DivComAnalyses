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

#' @title Plot taxa abundances over time
#' @author Sneha Sundar
#' @param physeq The phyloseq object
#' @param taxa_level The taxonomic level for which you want to plot abundances as a string (If you want to plot specific ASVs input "Species" as the argument)
#' @param taxa_to_plot The names of the taxa that you want to plot as a character vector 
#' @param time_column Name of the metadata column containing time information
#' @param other_columns_to_plot The names of metadata columns you want to additionally plot as a character vector(for example if you have total bacterial count for every sample)
#' @param axis_transform logical value indicating if you want to transform the axes prior to plotting 
#' @param transformation String indicating the transformation you want to perform (required if transform = TRUE); transformation chosen should be supported by the scales package; example of transformations: 'log10','log1p'
#' @param plot_total logical value indicating if you want to also see the total abundance of the taxa you selected in the plot (default: FALSE)
#' @param data_facet1 Name of a metadata column; category names will rows in the facet
#' @param data_facet2 Name of the metadata column; category names will columns in the facet
#' @param n_facet_col Number of facet rows
#' @param n_facet_row Number of facet columns
#' @note Make sure your phyloseq object is appropriately normalized (either to relative abundances or density normalization with qPCR data)
#' @return A list containing containing the dataframe used for plotting \code{plot.df} and the actual plot \code{plot}
#' @export
#' @examples
#' 
#'



plot_taxa_abundances_over_time<-function(physeq=physeq,
                                         taxa_level, 
                                         taxa_to_plot,
                                         time_column,
                                         other_columns_to_plot = NULL,
                                         axis_transform=TRUE,
                                         transformation = NULL,
                                         plot_total = FALSE,
                                         data_facet1=NULL,
                                         data_facet2=NULL,
                                         n_facet_row=NULL,
                                         n_facet_col=NULL){
  
  require(phyloseq)
  require(tidyverse)
  
  #if we want to plot specific ASVs
  if(taxa_level=="Species"){
    
    #make otu table filterable by ASV
    otu.table <- as(otu_table(physeq),"matrix") %>% data.frame(check.names=FALSE) %>% rownames_to_column("ASV")
    
    #add columns containing the percent of selected ASVs to sample data and make sure all the columns that should be numeric are numeric
    sample.data <- as(sample_data(physeq),"matrix") %>% data.frame(check.names=FALSE) %>% rownames_to_column("SampleID")
    
    #making the dataframe for plotting:
    
    #adding columns of the abundances of ASVs we want to plot
    for(i in taxa_to_plot){
      
      sample.data<-sample.data %>% mutate("{i}" :=  as.numeric(otu.table[otu.table$ASV==i,-1]))
    }
    #if `plot_total` == TRUE , add another column containing the total abundances of the ASVs we would like to plot
    if(plot_total ==TRUE){
      
      sample.data %>% mutate(ASV_total = 0) -> sample.data
      for(i in 1:length(taxa_to_plot)){
        sample.data %>% mutate(ASV_total = ASV_total+.data[[taxa_to_plot[i]]]) -> sample.data
      }
    }
    
    #make sure time column is numeric 
    sample.data %>% mutate("{time_column}":= as.numeric(get(time_column))) -> sample.data
    
    for(i in other_columns_to_plot){
      sample.data %>% mutate("{i}":= as.numeric(.data[[i]])) -> sample.data
    }
    
    #get the dataframe in long format for plotting with ggplot
    sample.data %>% 
      pivot_longer(cols = c(starts_with("ASV"),all_of(other_columns_to_plot)),names_to = "Data",values_to = "Abundances" ) %>% 
      arrange(get(time_column)) %>% 
      remove_missing(vars = "Abundances" ) -> plot.df
    
    #plot depending on whether we want a transformed axis or not
    if(axis_transform==TRUE){
      
      plot<-ggplot(data=plot.df,mapping=aes(color=Data,shape=Abundances==0))+
        geom_point(mapping=aes(x=.data[[time_column]],y=Abundances),alpha=0.8)+
        geom_path(mapping=aes(x=.data[[time_column]],y=Abundances,group=Data),alpha=0.8) +
        scale_y_continuous(trans=transformation)
    }else{
      plot<-ggplot(data=plot.df,mapping=aes(color=Data,shape=Abundances==0))+
        geom_point(mapping=aes(x=get(time_column),y=Abundances),alpha=0.8)+
        geom_path(mapping=aes(x=get(time_column),y=Abundances,group=Data),alpha=0.8) +
        geom_path(mapping=aes(x=get(time_column),y=Abundances,group=Data),alpha=0.8) 
      
    }
    
    #facet according to user defined input
    if(is.null(data_facet1) & !is.null(data_facet2)){
      plot<-plot + facet_wrap(~get(data_facet2),ncol=n_facet_col,nrow=n_facet_row) 
    }else if(is.null(data_facet2) & !is.null(data_facet1)){
      plot<- plot + facet_wrap(~get(data_facet1),ncol=n_facet_col,nrow=n_facet_row)
    }else if(!is.null(data_facet1) & !is.null(data_facet2)){
      plot<-plot+facet_wrap(get(data_facet1)~get(data_facet2),ncol=n_facet_col,nrow=n_facet_row)
    }
  }
  
  if(taxa_level!="Species"){
    
    physeq %>% tax_glom(taxrank = taxa_level) -> physeq
    
    #make otu table and tax table filterable by ASV
    otu.table <- as(otu_table(physeq),"matrix") %>% data.frame(check.names=FALSE) %>% rownames_to_column("ASV")
    tax.table <- as(tax_table(physeq),"matrix") %>% data.frame(check.names=FALSE) %>% rownames_to_column("ASV")
    #add columns containing the percent of selected ASVs to sample data and make sure all the columns that should be numeric are numeric
    sample.data <- as(sample_data(physeq),"matrix") %>% data.frame(check.names=FALSE) %>% rownames_to_column("SampleID")
    
    for(i in 1:length(taxa_to_plot)){
      ASV_to_plot<-tax.table$ASV[tax.table[[taxa_level]]==taxa_to_plot[i]]
      sample.data<-sample.data %>% mutate("{taxa_to_plot[i]}" :=  as.numeric(otu.table[otu.table$ASV==ASV_to_plot,-1]))
    }
    
    #if `plot_total` == TRUE , add another column containing the total abundances of the taxa we would like to plot
    if(plot_total ==TRUE){
      
      sample.data %>% mutate(Taxa_total = 0) -> sample.data
      for(i in 1:length(taxa_to_plot)){
        sample.data %>% mutate(Taxa_total = Taxa_total+.data[[taxa_to_plot[i]]]) -> sample.data
      }
      other_columns_to_plot <- c('Taxa_total',other_columns_to_plot)
    }
    
    sample.data %>% mutate("{time_column}":= as.numeric(get(time_column))) -> sample.data
    
    for(i in other_columns_to_plot){
      sample.data %>% mutate("{i}":= as.numeric(.data[[i]])) -> sample.data
    }
    
    sample.data %>% 
      pivot_longer(cols = c(all_of(taxa_to_plot),all_of(other_columns_to_plot)),names_to = "Data",values_to = "Abundances" ) %>%
      arrange(get(time_column)) %>% 
      remove_missing(vars = "Abundances" ) -> plot.df
    
    
    #plot depending on whether we want a transformed axis or not
    if(axis_transform==TRUE){
      
      plot<-ggplot(data=plot.df,mapping=aes(color=Data,shape=Abundances==0))+
        geom_point(mapping=aes(x=.data[[time_column]],y=Abundances),alpha=0.8)+
        geom_path(mapping=aes(x=.data[[time_column]],y=Abundances,group=Data),alpha=0.8) +
        scale_y_continuous(trans=transformation)
    }else{
      plot<-ggplot(data=plot.df,mapping=aes(color=Data,shape=Abundances==0))+
        geom_point(mapping=aes(x=get(time_column),y=Abundances),alpha=0.8)+
        geom_path(mapping=aes(x=get(time_column),y=Abundances,group=Data),alpha=0.8) +
        geom_path(mapping=aes(x=get(time_column),y=Abundances,group=Data),alpha=0.8) 
      
    }
    
    #facet according to user defined input
    if(is.null(data_facet1) & !is.null(data_facet2)){
      plot<-plot + facet_wrap(~get(data_facet2),ncol=n_facet_col,nrow=n_facet_row) 
    }else if(is.null(data_facet2) & !is.null(data_facet1)){
      plot<- plot + facet_wrap(~get(data_facet1),ncol=n_facet_col,nrow=n_facet_row)
    }else if(!is.null(data_facet1) & !is.null(data_facet2)){
      plot<-plot+facet_wrap(get(data_facet1)~get(data_facet2),ncol=n_facet_col,nrow=n_facet_row)
    }
    
  }
  
  return(out <- list("plot" = plot,
                     "plot.df" = plot.df))
  
  #code checks: what happens if the time column is a date?
  
}