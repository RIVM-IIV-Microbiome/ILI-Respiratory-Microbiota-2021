#' Calculate Fold Difference in Taxa Abundance
#'
#' @name calculateTaxaFoldDifference
#'
#' @details Calculate Log10 fold difference in abundance of taxa between two-
#'          groups. This code is modified from original code
#'          \url{https://github.com/microbiome/microbiome/blob/fa2c0de2fbe000da87be3c185972ed7f0f626591/inst/extdata/check_foldchange.R} Get the prevalence of taxa in \code{phyloseq} objects
#'          along with taxonomic classification and prevalence.
#'
#' @param x A phyloseq object
#'
#' @param group Vector with specifying the groups to compare.
#'              Only two-group comparisons are supported.
#'
#' @param sort Sort the results by descending order of fold difference
#'
#' @param paired Paired comparison (Default: FALSE)
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' # reduce size for example
#' ps1 <- subset_samples(FuentesIliGutData, ILI != "L2")
#' ps1 <- prune_taxa(taxa_sums(ps1) > 0, ps1)
#'
#' taxa_fd <- calculateTaxaFoldDifference(ps1, group="ILI")
#' # check
#' taxa_fd
#'
#' @return A tibble with taxa ids, taxonomic information,
#'         two-group prevalence and fold change values.
#'
#' @author Original Author: Leo Lahti. Adapted by Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2021). Utilities for microbiome analytics.
#' \url{https://github.com/microsud/biomeUtils}
#'
#' @importFrom tibble rownames_to_column
#' @importFrom phyloseq sample_names prune_samples
#' @importFrom microbiome abundances meta
#' @importFrom dplyr left_join relocate rename arrange desc filter mutate
#'
#' @export
NULL

calculateTaxaFoldDifference <- function (x, group, sort = FALSE, paired = FALSE) {
  
  #global vars
  FoldDifference <- NULL
  # Pick the grouping variable from sample metadata
  g <- group
  if (length(g) == 1) {
    g <- sample_data(x)[[g]]
    if (!is.factor(g)) {
      warning(paste("Converting the grouping variable", group, "into a factor."))
      g <- as.factor(g)
    }
    g <- droplevels(g)
    if (!length(levels(g)) == 2) {
      stop(paste("calculateTaxaFoldChange is currently implemented only for
                 two-group comparisons. The selected variable", group, "has",
                 length(unique(g)), "levels: ", paste(unique(g), collapse = "/")))
    }
  }
  if (is(x) == "phyloseq") {
    tx.ab <- abundances(x)
   
    #tx.ab <- abundances(x, "clr")
  }
  # Calculate fold changes
  if (paired) {
    fc <- apply(tx.ab, 1, function (xi) {spl <- split(xi, g);
    mean(spl[[2]] - spl[[1]], na.rm = TRUE)
    })
  } else {
    fc <- apply(tx.ab, 1, function (xi) {spl <- split(xi, g);
    mean(spl[[2]], na.rm = TRUE) - mean(spl[[1]], na.rm = TRUE)
    #log10(mean(spl[[2]], na.rm = TRUE)) - log10(mean(spl[[1]], na.rm = TRUE))
    })
  }
  
  fcdf <- as.data.frame(fc) %>%
    dplyr::rename(FoldDifference = fc) %>%
    rownames_to_column("FeatureID")
  
  
  lev.g1 <- levels(g)[1]
  lev.g2 <- levels(g)[2]
  
  g1.sams <- getSampleTibble(x,
                             select_rows = sample_names(x),
                             select_cols = group) %>%
    filter(.data[[group]] %in% lev.g1)
  
  g2.sams <- getSampleTibble(x,
                             select_rows = sample_names(x),
                             select_cols = group) %>%
    filter(.data[[group]] %in% lev.g2)
  
  lab.prev1 <- paste("Prevalence.",lev.g1, sep = "")
  lab.prev2 <- paste("Prevalence.",lev.g2, sep = "")
  
  prev.tb.g1 <- getPrevalence(prune_samples(sample_names(x) %in% g1.sams$SampleID, x),
                              return_rank= rank_names(x),
                              return_taxa = taxa_names(x))
  colnames(prev.tb.g1)[2] <- lab.prev1
  
  prev.tb.g2 <- getPrevalence(prune_samples(sample_names(x) %in% g2.sams$SampleID, x),
                              return_rank= rank_names(x),
                              return_taxa = taxa_names(x))
  colnames(prev.tb.g2)[2] <- lab.prev2
  
  prev.tb <- prev.tb.g1 %>%
    left_join(prev.tb.g2)
  
  fcdf <- fcdf %>%
    left_join(prev.tb, by= c(FeatureID = "Taxa")) %>%
    mutate(Enriched = ifelse(FoldDifference > 0, lev.g2,
                             ifelse(FoldDifference < 0 , lev.g1, "NoChange"))) # %>%
    #dplyr::relocate(rank_names(x), before = "FoldDifference")
  
  if (sort) {
    fcdf <- fcdf %>%
      arrange(desc(FoldDifference))
  }
  
  return(fcdf)
  
  }

###########################################################################
#' Get tibble
#'
#' @name getTibble
#'
#' @details Convert different \code{phyloseq} slots into tibbles.
#'
#' \code{getAbundanceTibble} gets the otu_table in tibble format.
#'
#' \code{getTaxaTibble} gets the taxa_table in tibble format.
#'
#' \code{getSampleTibble} gets the sample_sata in tibble format.
#'
#' @param x \code{phyloseq} object
#'
#' @param column_id A character
#'
#' @param select_rows Rows to return in output.
#'
#' @param select_cols Columns to return in output
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' otu_tib <- getAbundanceTibble(FuentesIliGutData)
#' tax_tib <- getTaxaTibble(FuentesIliGutData)
#' meta_tib <- getSampleTibble(FuentesIliGutData)
#' @return A tibble
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2020). Utilities for microbiome analytics.
#' \url{https://github.com/microsud/biomeUtils}
#'
#'
NULL


#' @rdname getTibble
#' @aliases getAbundanceTibble
#' @importFrom tibble rownames_to_column
#' @importFrom microbiome abundances
#' @importFrom dplyr as_tibble
#' @importFrom magrittr %>%
#' @importFrom rlang .data sym
#' @export
getAbundanceTibble <- function(x,
                               column_id = "FeatureID",
                               select_rows = NULL,
                               select_cols = NULL) {
  
  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }
  
  #column_id_1 <- column_id
  sel_rows <- .select_rows_abundance(x,select_rows)
  sel_cols <- .select_cols_abundance(x,select_cols)
  
  cols <- c(column_id, sel_cols)
  id.col <- sym(column_id)
  tib_dat <- abundances(x) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    rownames_to_column(column_id) %>%
    as_tibble() %>%
    filter(.data[[column_id]] %in% sel_rows)
  #filter(.data[[var]] %in% sel_rows)
  
  tib_dat <- tib_dat[,cols]
  return(tib_dat)
}


#' @rdname getTibble
#' @aliases getTaxaTibble
#' @importFrom tibble rownames_to_column
#' @importFrom phyloseq tax_table
#' @importFrom dplyr as_tibble
#' @importFrom methods as
#' @importFrom magrittr %>%
#' @importFrom rlang .data sym
#'
#' @export
getTaxaTibble <- function(x,
                          column_id = "FeatureID",
                          select_rows = NULL,
                          select_cols = NULL) {
  
  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }
  
  #column_id_1 <- column_id
  sel_rows <- .select_rows_taxonomy(x,select_rows)
  sel_cols <- .select_cols_taxonomy(x,select_cols)
  
  cols <- c(column_id, sel_cols)
  tib_dat <- tax_table(x) %>%
    as("matrix") %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    rownames_to_column(column_id) %>%
    as_tibble()%>%
    filter(.data[[column_id]] %in% sel_rows)
  
  tib_dat <- tib_dat[,cols]
  
  return(tib_dat)
}


#' @rdname getTibble
#' @aliases getSampleTibble
#' @importFrom tibble rownames_to_column
#' @importFrom microbiome meta
#' @importFrom dplyr as_tibble
#' @importFrom magrittr %>%
#' @importFrom rlang .data sym
#'
#' @export
getSampleTibble <- function(x,
                            column_id = "SampleID",
                            select_rows = NULL,
                            select_cols = NULL) {
  
  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }
  
  #column_id_1 <- column_id
  sel_rows <- .select_rows_sample(x,select_rows)
  sel_cols <- .select_cols_sample(x,select_cols)
  
  cols <- c(column_id, sel_cols)
  
  tib_dat <- meta(x) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    rownames_to_column(column_id) %>%
    as_tibble()%>%
    filter(.data[[column_id]] %in% sel_rows)
  
  tib_dat <- tib_dat[,cols]
  
  return(tib_dat)
}


###################### For getSampleTibble ######################
#' @importFrom phyloseq sample_variables

.select_cols_sample <- function(x, select_cols) {
  
  if(is.null(select_cols) || any(is.na(select_cols))){
    
    return(sample_variables(x))
    
  } else if(any(select_cols %in% sample_variables(x))) {
    
    return(select_cols)
    
  } else if(!is.null(select_cols) && !any(is.na(select_cols)) && !any(select_cols %in% sample_variables(x))){
    stop("For getSampleTibble, please provide valid sample_variables(x) in 'select_cols' ")
  }
  
}

#' @importFrom phyloseq sample_names
.select_rows_sample <- function(x, select_rows) {
  
  if(is.null(select_rows) || any(is.na(select_rows))){
    
    return(sample_names(x))
    
  } else if(any(select_rows %in% sample_names(x))) {
    
    return(select_rows)
    
  } else if(!is.null(select_rows) && !any(is.na(select_rows)) && !any(select_rows %in% sample_names(x))){
    stop("For getSampleTibble, please provide valid sample_names(x) in 'select_rows' ")
  }
  
}

###################### For getTaxaTibble ######################
#' @importFrom phyloseq rank_names
.select_cols_taxonomy <- function(x, select_cols) {
  
  if(is.null(select_cols) || any(is.na(select_cols)) ){
    
    return(rank_names(x))
    
  } else if(any(select_cols %in% rank_names(x))) {
    
    return(select_cols)
    
  } else if(!is.null(select_cols) && !any(is.na(select_cols)) && !any(select_cols %in% rank_names(x))){
    stop("For getTaxaTibble, provide valid rank_names(x) in 'select_cols' ")
  }
  
}

#' @importFrom phyloseq taxa_names
.select_rows_taxonomy <- function(x, select_rows) {
  
  #)
  if(is.null(select_rows) || any(is.na(select_rows))){
    
    return(taxa_names(x))
    
  } else if(any(select_rows %in% taxa_names(x))) {
    
    return(select_rows)
    
  } else if(!is.null(select_rows) && !any(is.na(select_rows)) && !any(select_rows %in% taxa_names(x))){
    stop("For getTaxaTibble, please provide valid taxa_names(x) in 'select_rows' ")
  }
  
}



###################### For getAbundanceTibble ######################
#' @importFrom phyloseq sample_names
.select_cols_abundance <- function(x, select_cols) {
  
  #
  if(is.null(select_cols) || any(is.na(select_cols))){
    
    return(sample_names(x))
    
  } else if(any(select_cols %in% sample_names(x))) {
    
    return(select_cols)
    
  } else if(!is.null(select_cols) && !any(is.na(select_cols)) && !any(select_cols %in% sample_names(x))){
    stop("For getAbundanceTibble, please provide valid sample_names(x) in 'select_cols' ")
  }
  
}

#' @importFrom phyloseq taxa_names
.select_rows_abundance <- function(x, select_rows) {
  
  #)
  if(is.null(select_rows) || any(is.na(select_rows))){
    
    return(taxa_names(x))
    
  } else if(any(select_rows %in% taxa_names(x))) {
    
    return(select_rows)
    
  } else if(!is.null(select_rows) && !any(is.na(select_rows)) && !any(select_rows %in% taxa_names(x))){
    stop("For getAbundanceTibble, please provide valid taxa_names(x) in 'select_rows' ")
  }
  
}

##################################################################
#' Get Prevalence and Taxonomy
#'
#' @name getPrevalence
#'
#' @details Get the prevalence of taxa in \code{phyloseq} objects
#'          along with taxonomic classification.
#'
#' @param x A phyloseq object
#'
#' @param return_rank Specify which taxonomic ranks to include in output.
#'                    Must be a character vector \code{phyloseq::rank_names()}
#'
#' @param return_taxa A specific list of taxa for which the values should
#'                    be returned. This can be used if user is not interested
#'                    in all the taxa in input \code{phyloseq}. Default is NULL
#'                    which returns all taxa. This list must match rows names if
#'                    otu_table in phyloseq has taxa_are_rows=TRUE or columns
#'                    names if otu_table in phyloseq has taxa_are_rows=FALSE
#'
#' @param sort Logical. Sort by prevalence value from higher to lower (Default=TRUE)
#'
#' @param ... Option to pass microbiome::prevalence
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' prev_tib <- getPrevalence(FuentesIliGutData,
#'                           return_rank= c("Family", "Genus"),
#'                           return_taxa = c("ASV4", "ASV17" , "ASV85", "ASV83"),
#'                           sort=TRUE)
#' head(prev_tib)
#'
#' @return A tibble with prevalence and taxonomy
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' \itemize{
#' \item{}{A Salonen et al. (2012) The adult intestinal core microbiota
#'         is determined by  analysis depth and health status.
#'        \emph{Clinical Microbiology and Infection}, 18(S4):16 20.
#'        \url{https://doi.org/10.1111/j.1469-0691.2012.03855.x}
#' }
#' \item{}{Lahti L, Shetty S (2012-2019). microbiome R package.
#'         statistical aspects of the study of the microbiome.
#'         \emph{BioConductor}.
#'         \url{https://doi.org/doi:10.18129/B9.bioc.microbiome}
#' }
#' }
#'
#'
#' @importFrom microbiome prevalence
#' @importFrom tibble rownames_to_column
#' @importFrom phyloseq rank_names tax_table
#' @importFrom dplyr arrange left_join select rename filter desc
#' @importFrom rlang syms
#'
#' @export
NULL

getPrevalence <- function(x,
                          return_rank= rank_names(x),
                          return_taxa = NULL,
                          sort = TRUE, ...) {
  
  # global vars
  Taxa <- NULL
  
  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }
  
  rank_nm <- .check_ranks(x, return_rank)
  rank_nm <- c("Taxa", rank_nm)
  
  tax_tib <- getTaxaTibble(x,
                           column_id = "Taxa",
                           select_cols = rank_names(x),
                           select_rows = taxa_names(x)) %>%
    dplyr::select(!!!syms(rank_nm))
  
  prev_tbl <- prevalence(x, ...) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Taxa") %>%
    dplyr::rename(prevalence = ".") %>%
    dplyr::left_join(tax_tib)
  
  #colnames(prev_tbl) <- c("Taxa", "prevalence", rank_nm)
  
  taxa_rt <- .check_taxa(x, return_taxa)
  
  prev_tbl <- prev_tbl %>%
    filter(Taxa %in% !!syms(taxa_rt))
  
  if (sort){
    
    prev_tbl <- prev_tbl %>%
      arrange(desc(prevalence))
    return(prev_tbl)
    
  } else {
    
    return(prev_tbl)
    
  }
  
  
  
}


#' @importFrom phyloseq rank_names
.check_ranks <- function(x, return_rank){
  
  if(!is.null(return_rank) || !is.na(return_rank) || any(return_rank %in% rank_names(x))){
    
    return(return_rank)
    
  } else if(is.null(return_rank) || is.na(return_rank)) {
    
    return_rank <- rank_names(x)
    return(return_rank)
    
  } else if(!is.null(return_rank) && !is.na(return_rank) && !any(return_rank %in% rank_names(x))){
    stop("Please provide valid taxonomic rank names in 'return_rank' ")
  }
  
}

#' @importFrom phyloseq taxa_names
.check_taxa <- function(x, return_taxa){
  
  if(!is.null(return_taxa) || !any(is.na(return_taxa))  || any(return_taxa %in% taxa_names(x))){
    
    return(return_taxa)
    
  } else if(is.null(return_taxa) || any(is.na(return_taxa))) {
    
    return(taxa_names(x))
    
  } else if(!is.null(return_taxa) && !any(is.na(return_taxa)) && !any(return_taxa %in% taxa_names(x))){
    stop("Please provide valid names in 'return_taxa' ")
  }
  
}


###
removeZeroTaxa <- function(x){
  return(prune_taxa(taxa_sums(x) >0, x))
}

