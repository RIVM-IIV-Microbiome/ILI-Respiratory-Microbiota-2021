# Root tree 
# source:https://john-quensen.com/r/unifrac-and-tree-roots/
pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <-
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup) }


transform_compositional = function(ps) {
  otu_table(ps) = otu_table(
    get_taxa(ps) / rowSums(get_taxa(ps)),
    taxa_are_rows = FALSE
  )
  ps
}


removeZeroTaxa <- function(x){
  return(prune_taxa(taxa_sums(x) >0, x))
}

# Calculate phylogenetic diversity  
# param x phyloseq object 
# retunr phyloseq object with results of PD added to sample_data()
phyloDiversity <- function(x){
  
  require(picante)
  
  meta_tb <- meta(x)
  otu_tb <- as.data.frame(abundances(x))
  asv.tree <- x@phy_tree
  df.pd <- pd(t(otu_tb), 
              asv.tree,
              include.root=T) # t(ou_table) transposes the table for use in picante and the tre file comes from the first code chunck we used to read tree file (see making a phyloseq object section).
  
  meta_tb <- cbind(meta_tb,df.pd)
  sample_data(x) <- sample_data(meta_tb)
  #sample_data(x)$Phylogenetic_Diversity <- df.pd$PD
  return(x)
}

################################################################################################
# Beta dispersion wrapper
# @param x phyloseq object
# @param group variable in metadata
beta_disp <- function(x, group="condition_status", index="bray"){
  dist3a <- phyloseq::distance(x, index)
  bet_di_anova3a <- anova(betadisper(dist3a, meta(x)[,group]))
  bet_di_anova3a
}

#################################################

get_taxa_tibble <- function(x){
  tax_table(x) %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("FeatureID")
}


################################################################################################
# Make pairs for stat_compare_means  
# Handy script for make comparison paris to use with ggpubr::stat_compare_means()
# @param x list of vectors 
# 
# example
# pir <- c("sad", "bad", "sad", "gad", "gad", "bad","sad", "bad", "sad", "gad", "gad", "bad")
# pl <- make_pairs(pir)
# pl

make_pairs <- function(x) {
  if (is.character(x)==TRUE) {
    #message("is char")
    var.lev <- unique(x)
  } else if (is.factor(x) == TRUE) {
    #message("is fac")
    var.lev <- levels(x)
  }
  # make a pairwise list that we want to compare.
  lev.pairs <- combn(seq_along(var.lev), 2, simplify = FALSE, FUN = function(i) var.lev[i] )
  return(lev.pairs)
}

################################################################################################
# Get a taxa table at specified level 
# @param x phyloseq object
# @param level.tax taxonomic level
get_tax_table <- function(x, level.tax, path.dir="ms/tabs/"){
  ps.lev <- aggregate_taxa(x, level = level.tax)
  phy_tab <- pseq_long(ps.lev, transform.counts = NULL) %>% 
    group_by(condition_status, unique) %>% 
    summarise(mean_abun = mean(abundance),
              sd_abun = sd(abundance))
  file.loc <- paste0(path.dir,level.tax,".txt")
  write.table(phy_tab, file.loc,sep="\t")
  return(phy_tab)
}

## Pairwise GUniFrac
calculate_gunifrac_anosim <- function(ps, dist, comparisons){
  
  ps.x <- ps
  ps.ls <- NULL
  in.dist <- dist
  comps_grp <- comparisons
  for(i in 1:length(comps_grp)){
    
    comps2 <- comps_grp[[i]]
    
    #print(comps2)
    group <- paste0(comps2[[1]], " vs ", comps2[[2]])
    print(group)
    
    sample_sel <- rownames(subset(meta(ps.x), condition_status %in% comps2))
    #sample_sel[1:10]
    ps.com <- prune_samples(sample_sel, ps.x)
    
    # ps.com <- subset_samples(ps.x, condition_status %in% comps_grp[[i]])
    
    ps.com <- prune_taxa((taxa_sums(ps.com) >0), ps.com)
    
    grps.compare <- sample_data(ps.com)$condition_status
    
    dist.mat <- in.dist[sample_names(ps.com),sample_names(ps.com)]
    
    anosim_unifrac <- vegan::anosim(as.dist(dist.mat), 
                                    grouping = grps.compare, 
                                    permutations=999)
    
    anosim.signif <-round(anosim_unifrac$signif, 5)
    anosim.R <-round(anosim_unifrac$statistic, 5)
    ps.ls <- rbind(ps.ls, c(group, anosim.R, anosim.signif))
  }
  ps.ls <- as.data.frame(ps.ls)
  
  colnames(ps.ls) <- c("Comparison","R", "Pval")
  return(ps.ls)
}


################################################################################################

plot_ordination_gg <- function(ps, ord){
  
  ordip <- plot_ordination(ps, ord, justDF = T)
  # Get axis 1 and 2 variation
  evals1 <- round(ord$values$Eigenvalues[1] / sum(ord$values$Eigenvalues) * 100, 2)
  evals2 <- round(ord$values$Eigenvalues[2] / sum(ord$values$Eigenvalues) * 100, 2)
  
  plot_ord_1 <- ggplot(ordip, aes(x = Axis.1, y = Axis.2,color = as.factor(condition_status))) + 
    # add layers
    # add lines
    geom_vline(xintercept = 0, alpha=0.5) + 
    geom_hline(yintercept = 0, alpha=0.5) +
    # add points
    geom_point(aes(color = as.factor(condition_status)), 
               size = 1, alpha=0.5) +
    # add gradient colors
    scale_color_manual("",values = dis_col) + 
    # add x and y labels
    xlab(paste("PCoA 1 (", evals1, "%)", sep = "")) +
    ylab(paste("PCoA 2 (", evals2, "%)", sep = "")) +
    theme_biome_utils() #labs(title = "ASV level") 
}


################################################################################################

################################################################################################

plot_centroids <- function(ps, ord){
  
  ordip <- plot_ordination(ps, ord, justDF = T)
  # Get axis 1 and 2 variation
  evals1 <- round(ord$values$Eigenvalues[1] / sum(ord$values$Eigenvalues) * 100, 2)
  evals2 <- round(ord$values$Eigenvalues[2] / sum(ord$values$Eigenvalues) * 100, 2)
  
  cents <- ordip %>% 
    group_by(condition_status) %>% 
    summarise(Axis.1a = mean(Axis.1),
              Axis.2a = mean(Axis.2),
              ster.1 = std.erF(Axis.1),
              ster.2 = std.erF(Axis.2))
  
  plot_ord_1 <- ggplot(cents, 
                       aes(x = Axis.1a, y = Axis.2a,
                           color = as.factor(condition_status))) + 
    geom_vline(xintercept = 0, alpha=0.5) + 
    geom_hline(yintercept = 0, alpha=0.5) +
    # add points
    geom_point(data=cents,aes(color = as.factor(condition_status)), 
               size = 3, alpha=0.5) +
    geom_errorbarh(data=cents,aes(xmin=Axis.1a-ster.1,xmax=Axis.1a+ster.1),width=0.1) +
    geom_errorbar(data=cents,aes(ymin=Axis.2a-ster.2,ymax=Axis.2a+ster.2),height=0.1) +
    #geom_errorbar(data=cents,aes(ymin=Axis.2a-ster.2,ymax=Axis.2a+ster.2),height=0.1) +
    scale_color_manual("",values = dis_col) + 
    # add x and y labels
    xlab(paste("PCoA 1 (", evals1, "%)", sep = "")) +
    ylab(paste("PCoA 2 (", evals2, "%)", sep = "")) +
    theme_biome_utils() #labs(title = "ASV level") 
}


################################################################################################

################################################################################################
# Convert phyloseq to long data frame  
# alternative to psmelt 
# @param pseq phyloseq object  
# @param transform.counts see microbiome::transform()
# example
# pl <- pseq_long(pseq, transform.counts = "NULL")
# head(p1)
pseq_long <- function(x, transform.counts = "NULL") {
  asv_tab <- asv_tab_lng <- tax_tab <- NULL
  if (is.null(transform.counts)) {
    x <- x
  }
  else if (transform.counts == "log10") {
    x <- transform(x, "log10")
  }
  else if (transform.counts == "Z-OTU") {
    x <- transform(x, "Z", "OTU")
  }
  else if (transform.counts == "Z-Sample") {
    x <- transform(x, "Z", "Sample")
  }
  else if (transform.counts == "compositional") {
    x <- transform(x, "compositional", "OTU")
  }
  else {
    stop("Please provide appropriate transformation")
  }
  asv_tab <- as.data.frame(abundances(x)) # get asvs/otus
  asv_tab$asv_id <- rownames(asv_tab) # add a new column for ids
  #head(asv_tab)
  asv_tab_lng <- asv_tab %>%
    reshape2::melt()
  #head(asv_tab_lng)
  colnames(asv_tab_lng) <- c("asv_id", "sample", "abundance")
  # tax_tab <- as.data.frame(tax_table(x)) # get taxonomy note: can be slow
  tax_tab <- as.data.frame(as(x@tax_table, "matrix")) # get taxonomy note: can be slow
  # tax_tab <- as.data.frame(tax_tab)
  tax_tab$asv_id <- rownames(tax_tab) # add a new column for ids
  asv_tax_tab <- asv_tab_lng %>%
    right_join(tax_tab, by = "asv_id")
  
  meta_tg <- meta(x)
  meta_tg$sample <- rownames(meta_tg)
  asv_tax_meta_tab <- asv_tax_tab %>%
    right_join(meta_tg, by = "sample")
  
  return(asv_tax_meta_tab)
}

##########################################################################################################3
# Calculate stability between two timepoints (vists) 
# @param ps_obj phyloseq object
# @param dist.method "bray", "canberra"
# This function is hardcoded requires specific variable names in metadata
stability <- function(ps_obj, dist.method="bray"){
  physeq <- NULL
  physeq <- ps_obj
  sample_data(physeq)$subject <- as.factor(sample_data(physeq)$participant_id)
  smp_df <- sample_data(physeq) %>%
    data.frame() %>%
    mutate_at(
      .funs = as.character, 
      .vars = c("Sample_Name", "participant_id", "SampleID")
    ) 
  
  #X <- t(abundances(physeq)) 
  #X <- asinh(X)
  #physeq <- phyloseq::phyloseq(
  #  otu_table(X, taxa_are_rows = FALSE),
  #  sample_data(physeq))
  physeq <- prune_taxa(taxa_sums(physeq) > 0 , physeq)
  physeq <- microbiome::transform(physeq, "compositional")
  #Y <- t(abundances(physeq))
  D <- phyloseq::distance(physeq, dist.method)
  #D <- cor(Y,Y, method= "spearman",use= "na.or.complete")
  
  df_dist <- reshape2::melt(
    as.matrix(D), 
    varnames = c("S1", "S2"),
    value.name = "dist") %>%
    #filter(S1 != S2) %>% 
    mutate_if(is.factor, as.character) %>% 
    left_join(smp_df, by = c("S1" = "Sample_Name")) %>%
    left_join(smp_df, by = c("S2" = "Sample_Name"),suffix = c("_1", "_2")) %>% 
    filter(S1 != S2) %>% 
    filter(participant_id_1 ==participant_id_2) %>%
    mutate(SubjectID =participant_id_1) %>%
    distinct(SubjectID, .keep_all = TRUE) %>% 
    dplyr::select(-participant_id_1, -participant_id_1) %>% 
    mutate(stability = 1-dist)
  return(df_dist)
}

##########################################################################################################3
# Calculate taxa proportional variability  
# @param pseq phyloseq object  
# pseq should be relative abundances
taxa_pv <- function(pseq) {
  
  sp_tab_tax <- NULL
  # proportional variability index (PV)
  pv_index <- function (Z){
    n = length(Z)
    pairs = combn(Z,2)
    min_z = apply(pairs,2, min)
    max_z = apply(pairs,2, max)
    z = 1- (min_z/max_z)
    PV=2*sum(z)/(n*(n-1))
    return(PV)
  }
  
  sp_tab_tax <- t(abundances(pseq))
  k <- NULL # for each of the taxa, we will add a seperate constant k.
  pv_val_tax <- NULL
  pv_df_tax <- NULL
  for(j in 1:ncol(sp_tab_tax)) {
    k = (1/100)*mean(sp_tab_tax[,j]) # 1% of mean rel abundance for each taxa
    #cd_val <- cd_index(sp_tab[,j]+ k)
    pv_val_tax <- pv_index(sp_tab_tax[,j] + k)
    #cv_val <- mean(sp_tab[,j] + k)/sd(sp_tab[,j] + k)
    pv_df_tax=rbind(pv_df_tax,c(colnames(sp_tab_tax)[j],pv_val_tax))
  }
  
  pv_df_tax <- as.data.frame(pv_df_tax)
  colnames(pv_df_tax) <- c("Taxa", "PVIndex")
  pv_df_tax$PVIndex <- as.numeric(as.character(pv_df_tax$PVIndex))
  return(pv_df_tax)
}

##########################################################################################################3
# Calculate taxa proportional variability bootstrapped 
# @param pseq phyloseq object  
# @param lower.conf=0.025
# @param upper.conf=0.975
# @param bs.iter=999
# pseq should be relative abundances
taxa_pv_boot <- function(pseq, lower.conf=0.025, upper.conf=0.975, bs.iter=99){
  rand_sams <- ps.sub <- txvp <- sx <- cis_df <- sxi <- NULL
  s <- c()
  for (i in seq_len(bs.iter)) {
    size_80 <- round(0.8*nsamples(pseq))
    rand_sams <- sample(sample_names(pseq), size= size_80, replace = TRUE)
    ps.sub <- prune_samples(sample_names(pseq) %in% rand_sams, pseq)
    #ps.sub <- prune_taxa(taxa_sums(ps.sub) > 0, ps.sub)
    txvp <- taxa_pv(ps.sub)
    #rownames(txvp) <- txvp$Taxa
    s[[i]] <- txvp
  }
  sx <- dplyr::bind_cols(s)
  sx <- sx %>%
    mutate(PVTax= Taxa) %>% 
    dplyr::select(starts_with("PV")) %>% 
    as.data.frame()
  rownames(sx) <- sx$PVTax
  sx$PVTax <- NULL
  
  cis <- c()
  for(tax in rownames(sx)){
    taxsp_lc <- quantile(sx[tax,],lower.conf,na.rm=TRUE)
    taxsp_uc <- quantile(sx[tax,],upper.conf,na.rm=TRUE)
    cis <- rbind(cis, c(tax,taxsp_lc, taxsp_uc))
  }
  cis_df <- as.data.frame(cis)
  colnames(cis_df) <- c("Taxa", "LowerCI", "UpperCI")
  sx$meanPV <- rowMeans(sx)
  #sx$Taxa <- rownames(sx)
  sxi <- cbind(sx,cis_df)
  return(sxi)
}


#' @title Heatmap using \code{\link{phyloseq-class}} and \code{\link{pheatmap}}
#' @description Plot heatmap using \code{\link{phyloseq-class}} object as input.
#' @param x \code{\link{phyloseq-class}} object.
#' @param subset.top either NA or number of Top OTUs to use for plotting.
#' @param transformation either 'log10', 'clr','Z', 'compositional', or NA
#' @param VariableA main variable of Interest.
#' @param heatcolors is the option for colors in \code{\link{pheatmap}}. Default is to use Spectral
#' @param ... Arguments to be passed \code{\link{pheatmap}}.
#' @return A \code{\link{pheatmap}} plot object.
#' @export
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @author Sudarshan A. Shetty (sudarshanshetty9@gmail.com)
#' @examples
#'
#' library(microbiomeutilities)
#' library(viridis)
#' library(RColorBrewer)
#' data("zackular2014")
#' ps0 <- zackular2014
#'
#' heat.sample <- plot_taxa_heatmap(ps0,
#'   subset.top = 20,
#'   VariableA = "DiseaseState",
#'   heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
#'   transformation = "log10"
#' )
#' @keywords visualization
plot_taxa_heatmap_2 <- function(x, subset.top, 
                                transformation,
                                VariableA, 
                                heatcolors = NULL, ...) {
  topOTU <- phyobj1 <- phyobj2 <- otu.mat <- NULL
  meta.tab <- select.meta <- color.heatmap <- NULL
  
  
  #x <- ps0
  
  x <- suppressWarnings(suppressMessages(format_to_besthit(x)))
  if (!is.null(subset.top)) {
    message(paste0("Top ", subset.top, " OTUs selected",
                   sep = " "
    ))
    topOTU <- top_taxa(x, n = subset.top)
  } else {
    stop("specify a number/value for subset.top")
  }
  
  prev.tx <- prevalence(x)[topOTU]
  
  if (transformation == "log10") {
    message("log10, if zeros in data then log10(1+x) will be used")
    phyobj1 <- prune_taxa(topOTU, x)
    message("First top taxa were selected and \nthen abundances tranformed to log10(1+X)")
    
    phyobj2 <- microbiome::transform(phyobj1, "log10")
  } else if (transformation == "compositional") {
    phyobjx <-  microbiome::transform(x, "compositional")
    phyobj2 <- prune_taxa(topOTU, phyobjx)
    
    #prev.tx <- prevalence(phyobj2)
    
    message("First converted to compositional \n then top taxa were selected")
  } else if (transformation == "Z-OTU") {
    phyobj1 <- prune_taxa(topOTU, x)
    
    #prev.tx <- prevalence(phyobj1)
    
    phyobj2 <-  microbiome::transform(phyobj1, "Z")
    message("First top taxa were selected and \nthen abundances tranformed to Z values")
  } else if (transformation == "clr") {
    phyobj1 <- prune_taxa(topOTU, x)
    
    #prev.tx <- prevalence(phyobj1)
    
    phyobj2 <-  microbiome::transform(phyobj1, "clr")
    message("First top taxa were selected and \nthen abundances tranformed to clr")
  } else if (!is.null(transformation)) {
    stop("specify a number for transformation, log10, compositional, Z-OTU, clr")
  }
  
  
  # format the taxonomy to incluse unique names
  # phyobj2 <- format_phyloseq(phyobj2)
  
  
  otu.mat <- abundances(phyobj2)
  meta.tab <- meta(phyobj2)
  
  # choose which variables of interest to include in
  # the heatmap
  select.meta <- subset(meta.tab, select = c(VariableA))
  
  
  if (is.null(heatcolors)) {
    color.heatmap <- brewer.pal(6, "Spectral")
  } else {
    color.heatmap <- heatcolors
  }
  
  newnames <- NULL
  newnames <- lapply(
    rownames(otu.mat),
    function(x) bquote(italic(.(x))))
  #row_df <- NULL
  #row_df <- as.data.frame(round(prev.tx*100, 2))
  #colnames(row_df) <- c("Prevalence")
  
  
  heatmap <- pheatmap::pheatmap(otu.mat,
                                labels_row = as.expression(newnames),
                                annotation_col = select.meta,
                                #annotation_colors = annotation_colors,
                                #annotation_row = row_df,
                                color = color.heatmap, ...)
  
  
  return(list("plot"=heatmap, "tax_tab"=otu.mat, "prev" = prev.tx))
  
}

