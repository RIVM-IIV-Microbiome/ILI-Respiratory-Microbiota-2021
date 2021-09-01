#' @title Create variable names
#' @param data.set dataset consisting of variables
#' @import graphics
#' @importFrom stats aggregate cor na.omit
#' @import utils
#' @import coin
#' @import dplyr
#' @import ggplot2
#' @import microbiome
#' @import phyloseq
#' @export

create.names <- function(data.set) {
  names.of.variables <- NULL
  for (j in (1:ncol(data.set))) {
    names.of.variables.j <- paste(unlist(strsplit(x = names(data.set)[j], split = "\\.")), collapse = " ")
    names.of.variables <- c(names.of.variables, paste(names.of.variables.j, sep = ""))
  }
  return(names.of.variables)
}

############################################################################################################

#' @title show.missing.by.variable
#' @param data.set dataset consisting of variables
#' @param plot Plot
#' @export

show.missing.by.variable <- function(data.set, plot) {
  missing.by.variable <- rep(0, ncol(data.set))
  names(missing.by.variable) <- names(data.set)
  for (j in (1:length(missing.by.variable))) {
    missing.by.variable[j] <- sum(is.na(data.set[, j]))
  }
  missing.by.variable <- missing.by.variable / nrow(data.set)
  if (plot) {
    par(mar = c(5, 12, 1, 2))
    barplot(missing.by.variable,
            horiz = TRUE, las = 1, xlab = "Proportion of missing observations",
            col = "lavender", cex.axis = 1.25, cex.lab = 1.25
    )
  }
  return(missing.by.variable)
}


############################################################################################################
#' @title check.strata
#' @param data.to.test dataset variable to specifically test
#' @param s Stratum that is created to check for in the analysis
#' @export
check.strata <- function(s, data.to.test) {
  aux.data.to.test <- subset(data.to.test, stratum == s)
  check <- min(length(unique(aux.data.to.test[, 1])), length(unique(aux.data.to.test[, 2])))
  return(check)
}
#

############################################################################################################

#' @title Wilson.interval
#' @param frequency frequency
#' @param n Number must be numeric
#' @param confidence confidence value
#' @export
#
Wilson.interval <- function(frequency,n,confidence){
  kappa <- qnorm(1-(1-confidence)/2)
  estimate <- frequency/n
  aux1 <- (frequency+(kappa^2)/2+kappa*sqrt(n)*
             sqrt(estimate*(1-estimate)+(kappa^2)/(4*n)))/(n+(kappa^2))
  aux2 <- (frequency+(kappa^2)/2-kappa*sqrt(n)*sqrt(estimate*(1-estimate)+(kappa^2)/(4*n)))/
    (n+(kappa^2))
  return(c(max(0.0,aux2),min(1.0,aux1)))
}
#
############################################################################################################



#' @title Transform to Fractions
#' @details Prepares the taxa abundance table for testing.
#' @param x \code{\link{phyloseq-class}} object
#' @param filter.abund remove taxa with less than this abundance in total dataset
#' @export

transform_to_fractions <- function(x,
                                   det.thres = 0.001,
                                   prev.thres = 5 / 100) {
  
  taxa_abund <- taxa_abund_tib <- sample_ldf_sub <- NULL
  Shannon <- value.reabund <- value <- variable <- SD <- NULL
  prop.positive <- CV <- entropy <- average <- NULL
  taxa_wide <- core_taxa <- taxa_ldf <- taxa_wide_core <- NULL
  ID <- sum.i <- NULL
  
  # x <- ps
  
  taxa_abund <- t(as.data.frame(microbiome::abundances(x)))
  taxa_abund_tib <- dplyr::as_tibble(taxa_abund, rownames = "ID")
  
  # names(asv_tab)[1] <- "ID"
  # get column of bacterial abundances
  # message("calculating ....")
  Shannon <- function(p) {
    entropy <- ifelse(p > 0, -p * log(p), 0)
    return(entropy)
  }
  
  taxa_ldf <- reshape2::melt(taxa_abund_tib)
  taxa_ldf <- taxa_ldf %>%
    group_by(ID) %>%
    mutate(
      average = mean(value, na.rm = TRUE),
      SD = sd(value, na.rm = TRUE),
      prop.positive = mean(value > 0, na.rm = TRUE),
      sum.i = sum(value, na.rm = TRUE),
      value.reabund = value / sum.i,
      entropy = sum(Shannon(value.reabund)),
      CV = SD / average
    ) %>%
    ungroup()
  
  # head(subset(asv_tab2, ID=="052_056_S56"))
  
  # other <- c("average","SD","prop.positive","entropy","CV")
  sample_ldf_sub <- taxa_ldf %>%
    dplyr::select(ID, average, SD, prop.positive, entropy, CV) %>%
    distinct(ID, .keep_all = T)
  
  taxa_wide <- taxa_ldf %>%
    tidyr::pivot_wider(
      id_cols = ID,
      names_from = variable,
      values_from = value.reabund
    ) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("ID")
  
  # Filter our data
  core_taxa <- core_members(t(taxa_wide),
                            detection = det.thres,
                            prevalence = prev.thres
  )
  taxa_wide_core <- taxa_wide[, core_taxa]
  
  # Join the selected taxa with properties of sample.
  taxa_wide_core <- taxa_wide_core %>%
    tibble::rownames_to_column("ID") %>%
    left_join(sample_ldf_sub, by = "ID")
  
  
  # message("done ....")
  
  return(taxa_wide_core)
}
##############################################################################3
#' @title Prepare Metadata
#' @details Prepares Metadata for Testing.
#' @param x \code{\link{phyloseq-class}} object
#' @param filter.abund remove taxa with less than this abundance in total dataset
#' @param sam.ids columns with sample ids matching rownames of otu_table(ps)
#' @param ref.group e.g. "no_ili" this will be coded as zero.
#' The ref group witin compare vars of interest
#' @param compare The column name with case-control description
#' @param confounder.vars confounders to create stratum e.g. c("gender", "age_yrs_oct2014")
#' @param var.interest variable of interest to test
#' @examples
#' # comord.vars <- c("BMI_2014","respiratory_disease_2014", "Roken_2014")
#' # med.vars <- c("ace_inhibitors_2014","arbs_2014","flu_shot_2014_2015")
#' # patho.vars <- c("influenza_2014_any","haemophilus_2014_any","coronavirus_2014")
#' # all.vars <- c(comord.vars,med.vars,patho.vars)
#' # phenotypic.data.1 <- prep_metadata_2(ps.genus.ctrl.ili,
#' # sam.ids = "seq_sam_name",
#' # ref.group="no_ili",
#' # compare= "condition_status",
#' # confounder.vars = c("gender", "age_yrs_oct2014"),
#' # var.interest=all.vars)
#'
#' @export

prep_metadata <- function(x,
                          sam.ids = NULL,
                          ref.group= NULL,
                          compare= NULL,
                          confounder.vars = NULL,
                          var.interest=NULL) {
  
  if(is.null(ref.group) | is.null(confounder.vars) | is.null(compare)){
    stop("Please supply main.group &  var.interest")
  }
  
  phenotypic.data <- NULL
  
  if(!is.null(var.interest)){
    var_ast_select <- c(sam.ids, compare, confounder.vars, var.interest)
  } else {
    var_ast_select <- c(sam.ids, compare, confounder.vars)
  }
  #x <-
  phenotypic.data <- microbiome::meta(x)[,var_ast_select]
  rownames(phenotypic.data) <- NULL
  names(phenotypic.data)[1] <- "ID"
  names(phenotypic.data)[2] <- "ILI"
  phenotypic.data <- base::transform(phenotypic.data,ILI=ifelse(ILI==ref.group,0,1))
  
  # covert all to numeric
  phenotypic.data <- phenotypic.data %>%
    mutate_if(is.character, as.numeric)
  # this needs consideration depending on data
  phenotypic.data <- phenotypic.data %>% replace(., is.na(.), -1)
  
  names(phenotypic.data)[4] <- "age"
  #phenotypic.data <- phenotypic.data[,-5] # remove age_group
  
  #range(phenotypic.data$age)
  age.breaks <- seq(60,90,5)
  phenotypic.data <- base::transform(phenotypic.data,
                                     age.group=cut(age,breaks=age.breaks,
                                                   right=TRUE,
                                                   include.lowest=TRUE))
  
  phenotypic.data$age.group <- as.character(phenotypic.data$age.group)
  #table(phenotypic.data$age.group)
  # Create strata of "confounders":
  confounders <- c("gender","age.group")
  name.of.stratification <- "Sex and age-group"
  #
  aux.data.set <- subset(phenotypic.data,select=confounders)
  aux.data.set <- base::transform(aux.data.set,stratum=rep(NA,nrow(aux.data.set)))
  for(i in (1:nrow(aux.data.set))){
    aux.data.set[i,]$stratum <- paste(aux.data.set[i,-length(aux.data.set)],collapse="/")
  }
  aux.data.set$stratum <- as.character(aux.data.set$stratum)
  # str(aux.data.set); head(aux.data.set)
  #
  phenotypic.data <- base::transform(phenotypic.data,stratum=aux.data.set$stratum)
  phenotypic.data$stratum <- as.character(phenotypic.data$stratum)
  #table(phenotypic.data$stratum)
  
  return(phenotypic.data)
}

##############################################################################3

#' @title Illustrate.ONE.association.BY.STRATUM
#' @details Tests based on the whole data set blocked by stratum.
#' This version uses the asymptotic distribution
#' @param columns dataset variable to specifically test
#' @param data.for.testing dataset to use for testing
#' @param TYPES variable properties, ordinal, categorical, continous, etc.
#' @param variables.of.interest variables to test
#' @export

test.ONE.association.0 <- function(columns, data.for.testing, TYPES, variables.of.interest) {
  # 	print(columns)
  # 	columns <- c(1,193); B <- 100000; Spearman.rather.than.AD <- TRUE; TYPES[columns]
  p.value <- test <- sample.characteristics <- Sign <- NA
  
  response.by.treatment <- na.omit(subset(data.for.testing,
                                          select = c(variables.of.interest[columns], "stratum")
  ))
  # 	head(response.by.treatment); str(response.by.treatment)
  
  checks.by.stratum <- t(sapply(response.by.treatment$stratum, check.strata, response.by.treatment))
  good.strata <- response.by.treatment$stratum[checks.by.stratum > 1]
  
  tested <- FALSE
  
  if (length(good.strata) > 1) {
    response.by.treatment <- subset(response.by.treatment, is.element(stratum, good.strata))
    names(response.by.treatment) <- c("treatment", "response", "stratum")
    rownames(response.by.treatment) <- NULL
    
    type.1 <- TYPES[columns[1]]
    type.2 <- TYPES[columns[2]]
    if ((type.1 == "binary" & type.2 != "binary") |
        (type.1 == "categorical" & (type.2 != "categorical" & type.2 != "binary")) |
        (type.1 == "ordinal" & (type.2 != "ordinal" & type.2 != "categorical" & type.2 != "binary"))) {
      response.by.treatment <- response.by.treatment[, c(2, 1, 3)]
      
      type.2 <- TYPES[columns[1]]
      type.1 <- TYPES[columns[2]]
    }
    names(response.by.treatment) <- c("response", "treatment", "stratum")
    head(response.by.treatment)
    response.by.treatment$stratum <- as.factor(response.by.treatment$stratum)
    
    if ((is.element(type.2, c("ordinal")) &
         is.element(type.1, c("dirac.and.continuous", "continuous", "ordinal"))) |
        (is.element(type.2, c("binary")) &
         is.element(type.1, c("ordinal")))) {
      independence.test <- independence_test(response ~ treatment | stratum,
                                             data = response.by.treatment,
                                             distribution = "asymptotic", teststat = "scalar"
      )
      p.value <- as.numeric(pvalue(independence.test))
      Sign <- sign(statistic(spearman_test(response ~ treatment | stratum,
                                           data = response.by.treatment,
                                           distribution = "asymptotic"
      )))
      list.of.samples <- split(response.by.treatment$response, f = response.by.treatment$treatment)
      test <- "Sum statistic"
      k <- length(list.of.samples)
      if (length(names(list.of.samples)) <= 5) {
        sample.characteristics <- NULL
        for (j in (1:k)) {
          sample.characteristics <- c(sample.characteristics, length(list.of.samples[[j]]))
        }
        sample.characteristics <- paste(sample.characteristics, collapse = "/")
      } else {
        sample.characteristics <- "too many"
      }
      tested <- TRUE
    }
    
    if (type.1 == "continuous" & type.2 == "continuous") {
      Spearman.test <- spearman_test(response ~ treatment | stratum,
                                     data = response.by.treatment,
                                     distribution = "asymptotic"
      )
      p.value <- as.numeric(pvalue(Spearman.test))
      Sign <- sign(statistic(Spearman.test))
      test <- "Spearman"
      sample.characteristics <- NA
      tested <- TRUE
    }
    
    if (is.element(type.2, c("binary", "categorical")) &
        is.element(type.1, c("dirac.and.continuous", "continuous")) & (tested == FALSE)) {
      KW.test <- kruskal_test(response ~ as.factor(treatment) | stratum,
                              data = response.by.treatment,
                              distribution = "asymptotic"
      )
      p.value <- as.numeric(pvalue(KW.test))
      if (type.2 == "binary") {
        Sign <- sign(cor(
          y = response.by.treatment$response,
          x = response.by.treatment$treatment
        ))
      } else {
        Sign <- 0
      }
      list.of.samples <- split(response.by.treatment$response, f = response.by.treatment$treatment)
      test <- "Kruskal-Wallis"
      k <- length(list.of.samples)
      sample.characteristics <- NULL
      for (j in (1:k)) {
        sample.characteristics <- c(sample.characteristics, length(list.of.samples[[j]]))
      }
      sample.characteristics <- paste(sample.characteristics, collapse = "/")
      tested <- TRUE
    }
    
    if (is.element(type.1, c("binary", "categorical", "ordinal")) &
        is.element(type.2, c("binary", "categorical", "ordinal")) & (tested == FALSE)) {
      CHM.test <- cmh_test(as.factor(response) ~ as.factor(treatment) | stratum,
                           data = response.by.treatment, distribution = "asymptotic"
      )
      p.value <- as.numeric(pvalue(CHM.test))
      if ((type.1 == "binary" & type.2 != "categorical") | (type.2 == "binary" & type.1 != "categorical")) {
        Sign <- sign(cor(
          y = response.by.treatment$response,
          x = response.by.treatment$treatment
        ))
      } else {
        Sign <- 0
      }
      test <- "Cochran-Mantel-Haenszel"
      list.of.samples <- split(response.by.treatment$response, f = response.by.treatment$treatment)
      k <- length(list.of.samples)
      sample.characteristics <- NULL
      for (j in (1:k)) {
        sample.characteristics <- c(sample.characteristics, length(list.of.samples[[j]]))
      }
      sample.characteristics <- paste(sample.characteristics, collapse = "/")
      tested <- TRUE
    }
  }
  
  if (tested == FALSE) {
    print(paste("Did not test ", paste(columns, collapse = ","), "!", sep = ""))
  }
  
  return(c(
    names(data.for.testing)[columns], as.character(p.value), Sign, test,
    sample.characteristics
  ))
}
#


##############################################################################3
#' @title Test one association
#' @details Tests based on the whole data set blocked by stratum
#' @param columns dataset variable to specifically test
#' @param data.for.testing dataset to use for testing
#' @param TYPES variable properties, ordinal, categorical, continous, etc.
#' @param variables.of.interest list/vector of variables matching input colnames to test
#' @param B.i resample values (numeric)
#' @export

test.ONE.association <- function(columns, data.for.testing, TYPES, variables.of.interest, B.i) {
  # 	print(columns)
  # 	columns <- c(1,193); B <- 100000; Spearman.rather.than.AD <- TRUE; TYPES[columns]
  p.value <- test <- sample.characteristics <- Sign <- NA
  
  response.by.treatment <- na.omit(subset(data.for.testing,
                                          select = c(variables.of.interest[columns], "stratum")
  ))
  # 	head(response.by.treatment); str(response.by.treatment)
  
  checks.by.stratum <- t(sapply(response.by.treatment$stratum, check.strata, response.by.treatment))
  good.strata <- response.by.treatment$stratum[checks.by.stratum > 1]
  head(good.strata)
  tested <- FALSE
  
  if (length(good.strata) > 1) {
    response.by.treatment <- subset(response.by.treatment, is.element(stratum, good.strata))
    names(response.by.treatment) <- c("treatment", "response", "stratum")
    
    rownames(response.by.treatment) <- NULL
    
    type.1 <- TYPES[columns[1]]
    type.2 <- TYPES[columns[2]]
    if ((type.1 == "binary" & type.2 != "binary") |
        (type.1 == "categorical" & (type.2 != "categorical" & type.2 != "binary")) |
        (type.1 == "ordinal" & (type.2 != "ordinal" & type.2 != "categorical" & type.2 != "binary"))) {
      response.by.treatment <- response.by.treatment[, c(2, 1, 3)]
      type.2 <- TYPES[columns[1]]
      type.1 <- TYPES[columns[2]]
    }
    names(response.by.treatment) <- c("response", "treatment", "stratum")
    response.by.treatment$stratum <- as.factor(response.by.treatment$stratum)
    
    if ((is.element(type.2, c("ordinal")) &
         is.element(type.1, c("dirac.and.continuous", "continuous", "ordinal"))) |
        (is.element(type.2, c("binary")) &
         is.element(type.1, c("ordinal")))) {
      independence.test <- independence_test(response ~ treatment | stratum,
                                             data = response.by.treatment,
                                             distribution = approximate(nresample = B.i), teststat = "scalar"
      )
      p.value <- as.numeric(pvalue(independence.test))
      Sign <- sign(statistic(spearman_test(response ~ treatment | stratum,
                                           data = response.by.treatment,
                                           distribution = "asymptotic"
      )))
      list.of.samples <- split(response.by.treatment$response, f = response.by.treatment$treatment)
      test <- "Sum statistic"
      k <- length(list.of.samples)
      if (length(names(list.of.samples)) <= 5) {
        sample.characteristics <- NULL
        for (j in (1:k)) {
          sample.characteristics <- c(sample.characteristics, length(list.of.samples[[j]]))
        }
        sample.characteristics <- paste(sample.characteristics, collapse = "/")
      } else {
        sample.characteristics <- "too many"
      }
      tested <- TRUE
    }
    
    if (type.1 == "continuous" & type.2 == "continuous") {
      Spearman.test <- spearman_test(response ~ treatment | stratum,
                                     data = response.by.treatment,
                                     distribution = approximate(nresample = B.i)
      )
      p.value <- as.numeric(pvalue(Spearman.test))
      Sign <- sign(statistic(Spearman.test))
      test <- "Spearman"
      sample.characteristics <- NA
      tested <- TRUE
    }
    
    if (is.element(type.2, c("binary", "categorical")) &
        is.element(type.1, c("dirac.and.continuous", "continuous")) & (tested == FALSE)) {
      KW.test <- kruskal_test(response ~ as.factor(treatment) | stratum,
                              data = response.by.treatment,
                              distribution = approximate(nresample = B.i)
      )
      p.value <- as.numeric(pvalue(KW.test))
      if (type.2 == "binary") {
        Sign <- sign(cor(
          y = response.by.treatment$response,
          x = response.by.treatment$treatment
        ))
      } else {
        Sign <- 0
      }
      list.of.samples <- split(response.by.treatment$response, f = response.by.treatment$treatment)
      test <- "Kruskal-Wallis"
      k <- length(list.of.samples)
      sample.characteristics <- NULL
      for (j in (1:k)) {
        sample.characteristics <- c(sample.characteristics, length(list.of.samples[[j]]))
      }
      sample.characteristics <- paste(sample.characteristics, collapse = "/")
      tested <- TRUE
    }
    
    if (is.element(type.1, c("binary", "categorical", "ordinal")) &
        is.element(type.2, c("binary", "categorical", "ordinal")) & (tested == FALSE)) {
      CHM.test <- cmh_test(as.factor(response) ~ as.factor(treatment) | stratum,
                           data = response.by.treatment, distribution = approximate(nresample = B.i)
      )
      p.value <- as.numeric(pvalue(CHM.test))
      if ((type.1 == "binary" & type.2 != "categorical") | (type.2 == "binary" & type.1 != "categorical")) {
        Sign <- sign(cor(
          y = response.by.treatment$response,
          x = response.by.treatment$treatment
        ))
      } else {
        Sign <- 0
      }
      test <- "Cochran-Mantel-Haenszel"
      list.of.samples <- split(response.by.treatment$response, f = response.by.treatment$treatment)
      k <- length(list.of.samples)
      sample.characteristics <- NULL
      for (j in (1:k)) {
        sample.characteristics <- c(sample.characteristics, length(list.of.samples[[j]]))
      }
      sample.characteristics <- paste(sample.characteristics, collapse = "/")
      tested <- TRUE
    }
  }
  
  if (tested == FALSE) {
    print(paste("Did not test ", paste(columns, collapse = ","), "!", sep = ""))
  }
  
  return(c(
    names(data.for.testing)[columns], as.character(p.value), Sign, test,
    sample.characteristics
  ))
}
#




##############################################################################3

#' @title Illustrate.ONE.association
#' @details Plots association
#' @param columns dataset variable to specifically test
#' @param data.for.plotting dataset to use for plotting
#' @param types variable properties, ordinal, categorical, continous, etc.
#' @param log.scale use log-scale to plot abundances
#' @export

illustrate.ONE.association <- function(columns, data.for.plotting, types, log.scale) {
  # print(columns)
  # 	columns <- columns.r; data.for.plotting <- data.for.testing
  
  plotted <- FALSE
  colour.palette <- c(
    "lavender", "lightblue", "cornflowerblue", "aquamarine2",
    "lightgreen", "olivedrab2", "palegreen1", "dodgerblue", "cyan",
    "bisque2", "burlywood2", "darkgoldenrod3",
    "bisque4", "black", "blueviolet", "darkorchid4", "blue",
    "cyan2", "aquamarine4", "chartreuse3", "darkolivegreen2", "darkseagreen3", "grey",
    "darkkhaki", "darkorange4", "brown3", "red", "darksalmon", "orange1", "yellow", "yellow3",
    "springgreen2", "lavender", "sienna2", "pink3", "khaki1",
    "darkmagenta", "azure2"
  )
  
  data.to.plot <- na.omit(data.for.plotting[, columns])
  type.1 <- types[columns[1]]
  type.2 <- types[columns[2]]
  types[columns]
  if ((type.2 == "binary" & type.1 != "binary") | (type.2 == "categorical" & type.1 != "binary")) {
    data.to.plot <- data.to.plot[, c(2, 1)]
    type.2 <- types[columns[1]]
    type.1 <- types[columns[2]]
  }
  good.names <- create.names(data.to.plot)
  names(data.to.plot) <- c("x", "y")
  
  if (type.1 == "dirac.and.continuous") {
    if (sum(data.to.plot$x == 0) <= 1) {
      type.1 <- "continuous"
    }
  }
  if (type.2 == "dirac.and.continuous") {
    if (sum(data.to.plot$y == 0) <= 1) {
      type.2 <- "continuous"
    }
  }
  
  if (is.element(type.1, c("categorical", "binary")) & type.2 == "dirac.and.continuous") {
    par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    aux.table <- table(data.to.plot[, 2])
    atom <- as.numeric(dimnames(aux.table)[[1]][which.max(as.vector(aux.table))])
    aux.data.to.plot <- base::transform(data.to.plot,
                                        y = ifelse(y > atom, paste(">", atom, sep = ""),
                                                   ifelse(y < atom, paste("<", atom, sep = ""), paste("=", atom, sep = ""))
                                        )
    )
    contingency.table <- table(aux.data.to.plot)
    frequency.table <- contingency.table / rowSums(contingency.table)
    barplot(t(frequency.table),
            beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
            ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)),
            xlab = good.names[1]
    )
    legend(
      x = "top", legend = colnames(frequency.table), horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:nrow(frequency.table)], title = good.names[2]
    )
    aux.data.to.plot <- subset(data.to.plot, y != atom)
    if (length(unique(aux.data.to.plot$x)) == 1) {
      x.label <- paste(good.names[1], "=", unique(aux.data.to.plot$x), collapse = "")
    } else {
      x.label <- good.names[1]
    }
    if (log.scale) {
      boxplot(log(y) ~ x,
              data = aux.data.to.plot, col = "lavender", xlab = x.label,
              ylab = paste("log(", good.names[2], ")", collapse = "")
      )
      stripchart(log(y) ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(y ~ x,
              data = aux.data.to.plot, col = "lavender", xlab = x.label,
              ylab = good.names[2]
      )
      stripchart(y ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    }
    plotted <- TRUE
  }
  
  if (is.element(type.1, c("categorical", "binary")) &
      is.element(type.2, c("categorical", "binary"))) {
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    contingency.table <- table(data.to.plot)
    frequency.table <- contingency.table / rowSums(contingency.table)
    barplot(t(frequency.table),
            beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
            ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)),
            xlab = good.names[1]
    )
    legend(
      x = "top", legend = dimnames(frequency.table)$y, horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:ncol(frequency.table)], title = good.names[2]
    )
    plotted <- TRUE
  }
  
  if (is.element(type.1, c("categorical", "binary", "ordinal")) & type.2 == "continuous") {
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    if (log.scale) {
      boxplot(log(y) ~ x,
              data = data.to.plot, col = "lavender", xlab = good.names[1],
              ylab = paste("log(", good.names[2], ")", collapse = "")
      )
      stripchart(log(y) ~ x, data = data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(y ~ x,
              data = data.to.plot, col = "lavender", xlab = good.names[1],
              ylab = good.names[2]
      )
    }
    plotted <- TRUE
  }
  
  if (is.element(type.2, c("categorical", "binary", "ordinal")) & type.1 == "continuous") {
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    if (log.scale) {
      boxplot(log(x) ~ y,
              data = data.to.plot, col = "lavender", xlab = good.names[2],
              ylab = paste("log(", good.names[1], ")", collapse = "")
      )
      stripchart(log(x) ~ y, data = data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(x ~ y,
              data = data.to.plot, col = "lavender", xlab = good.names[2],
              ylab = good.names[1]
      )
      stripchart(x ~ y, data = data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    }
    plotted <- TRUE
  }
  
  if ((type.1 == "ordinal" & type.2 == "ordinal") |
      (type.1 == "binary" & type.2 == "ordinal") | (type.1 == "categorical" & type.2 == "ordinal")) {
    aux.table.1 <- table(data.to.plot[, 1])
    aux.table.2 <- table(data.to.plot[, 2])
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    if (length(aux.table.2) > length(aux.table.1)) {
      aux.means <- aggregate(y ~ x, data = data.to.plot, mean)
      boxplot(y ~ x, data = data.to.plot, col = "lavender", xlab = good.names[1], ylab = good.names[2])
      points(y ~ as.factor(x),
             data = aux.means, col = "red", xlab = good.names[1], ylab = good.names[2],
             cex = 1.25, pch = 19
      )
    } else {
      aux.means <- aggregate(x ~ y, data = data.to.plot, mean)
      boxplot(x ~ y, data = data.to.plot, col = "lavender", xlab = good.names[2], ylab = good.names[1])
      points(x ~ as.factor(y),
             data = aux.means, col = "red", pch = 19, cex = 1.25,
             xlab = good.names[2], ylab = good.names[1]
      )
    }
    plotted <- TRUE
  }
  
  if (type.1 == "continuous" & type.2 == "continuous") {
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    plot(y ~ x, data = data.to.plot, col = "blue", xlab = good.names[1], ylab = good.names[2])
    plotted <- TRUE
  }
  
  if (type.1 == "ordinal" & type.2 == "dirac.and.continuous") {
    par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    aux.table <- table(data.to.plot[, 2])
    atom <- as.numeric(dimnames(aux.table)[[1]][which.max(as.vector(aux.table))])
    aux.data.to.plot <- base::transform(data.to.plot,
                                        y = ifelse(y > atom, paste(">", atom, sep = ""),
                                                   ifelse(y < atom, paste("<", atom, sep = ""), paste("=", atom, sep = ""))
                                        )
    )
    contingency.table <- t(table(aux.data.to.plot))
    frequency.table <- contingency.table / rowSums(contingency.table)
    barplot(t(frequency.table),
            beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
            ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)),
            xlab = good.names[2]
    )
    legend(
      x = "top", legend = colnames(frequency.table), horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:ncol(frequency.table)], title = good.names[1]
    )
    aux.data.to.plot <- subset(data.to.plot, y != atom)
    if (log.scale) {
      boxplot(log(y) ~ x,
              data = aux.data.to.plot, col = "lavender", xlab = good.names[1],
              ylab = paste("log(", good.names[2], ")", collapse = "")
      )
      stripchart(log(y) ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(y ~ x,
              data = aux.data.to.plot, col = "lavender", xlab = good.names[1],
              ylab = good.names[2]
      )
      stripchart(y ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    }
    plotted <- TRUE
  }
  
  if (type.2 == "ordinal" & type.1 == "dirac.and.continuous") {
    par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    aux.table <- table(data.to.plot[, 1])
    atom <- as.numeric(dimnames(aux.table)[[1]][which.max(as.vector(aux.table))])
    aux.data.to.plot <- base::transform(data.to.plot,
                                        x = ifelse(x > atom, paste(">", atom, sep = ""),
                                                   ifelse(x < atom, paste("<", atom, sep = ""), paste("=", atom, sep = ""))
                                        )
    )
    contingency.table <- table(aux.data.to.plot)
    frequency.table <- contingency.table / rowSums(contingency.table)
    barplot(t(frequency.table),
            beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
            ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)),
            xlab = good.names[1]
    )
    legend(
      x = "top", legend = colnames(frequency.table), horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:ncol(frequency.table)], title = good.names[2]
    )
    aux.data.to.plot <- subset(data.to.plot, x != atom)
    if (log.scale) {
      boxplot(log(x) ~ y,
              data = aux.data.to.plot, col = "lavender", xlab = good.names[2],
              ylab = paste("log(", good.names[1], ")", collapse = "")
      )
      stripchart(log(x) ~ y, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(x ~ y,
              data = aux.data.to.plot, col = "lavender", xlab = good.names[2],
              ylab = good.names[1]
      )
      stripchart(x ~ y, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    }
    plotted <- TRUE
  }
  
  if (type.1 == "dirac.and.continuous" & type.2 == "dirac.and.continuous") {
    aux.table.1 <- table(data.to.plot[, 1])
    aux.table.2 <- table(data.to.plot[, 2])
    atom.1 <- as.numeric(dimnames(aux.table.1)[[1]][which.max(as.vector(aux.table.1))])
    atom.2 <- as.numeric(dimnames(aux.table.2)[[1]][which.max(as.vector(aux.table.2))])
    aux.data.to.plot <- base::transform(data.to.plot,
                                        x = ifelse(x > atom.1, paste(">", atom.1, sep = ""),
                                                   ifelse(x < atom.1, paste("<", atom.1, sep = ""), paste("=", atom.1, sep = ""))
                                        ),
                                        y = ifelse(y > atom.2, paste(">", atom.2, sep = ""),
                                                   ifelse(y < atom.2, paste("<", atom.2, sep = ""), paste("=", atom.2, sep = ""))
                                        )
    )
    contingency.table <- table(aux.data.to.plot)
    frequency.table <- contingency.table / rowSums(contingency.table)
    aux.data.to.plot <- subset(data.to.plot, x != atom.1 & y != atom.2)
    if (nrow(aux.data.to.plot) >= 1) {
      par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    }
    barplot(t(frequency.table),
            beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
            ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)),
            xlab = good.names[1]
    )
    legend(
      x = "top", legend = colnames(frequency.table), horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:ncol(frequency.table)], title = good.names[2]
    )
    if (nrow(aux.data.to.plot) >= 1) {
      plot(y ~ x, data = aux.data.to.plot, col = "blue", xlab = good.names[1], ylab = good.names[2])
      stripchart(y ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    }
    plotted <- TRUE
  }
  
  if (type.1 == "dirac.and.continuous" & type.2 == "continuous") {
    par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    aux.table <- table(data.to.plot[, 1])
    atom <- as.numeric(dimnames(aux.table)[[1]][which.max(as.vector(aux.table))])
    aux.data.to.plot <- base::transform(data.to.plot,
                                        x = ifelse(x > atom, paste(">", atom, sep = ""),
                                                   ifelse(x < atom, paste("<", atom, sep = ""), paste("=", atom, sep = ""))
                                        )
    )
    if (log.scale) {
      boxplot(log(y) ~ x,
              data = aux.data.to.plot, col = "lavender", xlab = good.names[1],
              ylab = paste("log(", good.names[2], ")", collapse = "")
      )
      stripchart(log(y) ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(y ~ x,
              data = aux.data.to.plot, col = "lavender", xlab = good.names[1],
              ylab = good.names[2]
      )
      stripchart(y ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
      aux.data.to.plot <- subset(data.to.plot, x != atom)
      plot(y ~ x, data = aux.data.to.plot, col = "blue", xlab = good.names[1], ylab = good.names[2])
    }
    plotted <- TRUE
  }
  
  if (type.2 == "dirac.and.continuous" & type.1 == "continuous") {
    par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    aux.table <- table(data.to.plot[, 2])
    atom <- as.numeric(dimnames(aux.table)[[1]][which.max(as.vector(aux.table))])
    aux.data.to.plot <- base::transform(data.to.plot,
                                        y = ifelse(y > atom, paste(">", atom, sep = ""),
                                                   ifelse(y < atom, paste("<", atom, sep = ""), paste("=", atom, sep = ""))
                                        )
    )
    if (log.scale) {
      boxplot(log(x) ~ y,
              data = aux.data.to.plot, col = "lavender", xlab = good.names[2],
              ylab = paste("log(", good.names[1], ")", collapse = "")
      )
      stripchart(log(x) ~ y, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
      aux.data.to.plot <- subset(data.to.plot, y != atom)
      plot(x ~ y, data = aux.data.to.plot, col = "blue", xlab = good.names[2], ylab = good.names[1])
    } else {
      boxplot(x ~ y,
              data = aux.data.to.plot, col = "lavender", xlab = good.names[2],
              ylab = good.names[1]
      )
      stripchart(x ~ y, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
      aux.data.to.plot <- subset(data.to.plot, y != atom)
      plot(x ~ y, data = aux.data.to.plot, col = "blue", xlab = good.names[2], ylab = good.names[1])
    }
    plotted <- TRUE
  }
  
  if (((type.1 == "binary" & type.2 == "ordinal") | (type.1 == "categorical" & type.2 == "ordinal")) &
      (plotted == FALSE)) {
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    contingency.table <- table(data.to.plot)
    frequency.table <- contingency.table / rowSums(contingency.table)
    barplot(t(frequency.table),
            beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
            ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)),
            xlab = good.names[1]
    )
    legend(
      x = "top", legend = dimnames(frequency.table)$y, horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:ncol(frequency.table)], title = good.names[2]
    )
    plotted <- TRUE
  }
  
  if (plotted == FALSE) {
    print(paste("Did not plot ", paste(columns, collapse = ","), "!", sep = ""))
  }
}



##############################################################################3
#' @title Illustrate.ONE.association.BY.STRATUM
#' @details Plots association blocked by stratum
#' @param columns dataset variable to specifically test
#' @param data.for.plotting dataset to use for plotting
#' @param TYPES variable properties, ordinal, categorical, continous, etc.
#' @param stratum stratum to check
#' @param x.limits x-axis limits for plotting
#' @param y.limits y-axis limits for plotting
#' @param log.scale use log-scale to plot abundances
#' @export


illustrate.ONE.association.BY.STRATUM <- function(data.for.plotting, TYPES, stratum, x.limits, y.limits,
                                                  log.scale) {
  columns <- c(1, 2)
  plotted <- FALSE
  colour.palette <- c(
    "lavender", "lightblue", "cornflowerblue", "aquamarine2",
    "lightgreen", "olivedrab2", "palegreen1", "dodgerblue", "cyan",
    "bisque2", "burlywood2", "darkgoldenrod3",
    "bisque4", "black", "blueviolet", "darkorchid4", "blue",
    "cyan2", "aquamarine4", "chartreuse3", "darkolivegreen2", "darkseagreen3", "grey",
    "darkkhaki", "darkorange4", "brown3", "red", "darksalmon", "orange1", "yellow", "yellow3",
    "springgreen2", "lavender", "sienna2", "pink3", "khaki1",
    "darkmagenta", "azure2"
  )
  
  data.to.plot <- na.omit(data.for.plotting[, columns])
  type.1 <- TYPES[columns[1]]
  type.2 <- TYPES[columns[2]]
  TYPES[columns]
  if ((type.2 == "binary" & type.1 != "binary") | (type.2 == "categorical" & type.1 != "binary")) {
    data.to.plot <- data.to.plot[, c(2, 1)]
    type.2 <- TYPES[columns[1]]
    type.1 <- TYPES[columns[2]]
    aux.x.limits <- x.limits
    x.limits <- y.limits
    y.limits <- aux.x.limits
  }
  good.names <- create.names(data.to.plot)
  names(data.to.plot) <- c("x", "y")
  
  if (type.1 == "dirac.and.continuous") {
    if (sum(data.to.plot$x == 0) <= 1) {
      type.1 <- "continuous"
    }
  }
  if (type.2 == "dirac.and.continuous") {
    if (sum(data.to.plot$y == 0) <= 1) {
      type.2 <- "continuous"
    }
  }
  
  if (is.element(type.1, c("categorical", "binary")) & type.2 == "dirac.and.continuous") {
    par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    aux.table <- table(data.to.plot[, 2])
    atom <- as.numeric(dimnames(aux.table)[[1]][which.max(as.vector(aux.table))])
    aux.data.to.plot <- base::transform(data.to.plot,
                                        y = ifelse(y > atom, paste(">", atom, sep = ""),
                                                   ifelse(y < atom, paste("<", atom, sep = ""), paste("=", atom, sep = ""))
                                        )
    )
    contingency.table <- table(aux.data.to.plot)
    frequency.table <- contingency.table / rowSums(contingency.table)
    barplot(t(frequency.table),
            beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
            ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)),
            xlab = good.names[1], main = stratum
    )
    legend(
      x = "top", legend = colnames(frequency.table), horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:nrow(frequency.table)], title = good.names[2]
    )
    aux.data.to.plot <- subset(data.to.plot, y != atom)
    if (length(unique(aux.data.to.plot$x)) == 1) {
      x.label <- paste(good.names[1], "=", unique(aux.data.to.plot$x), collapse = "")
    } else {
      x.label <- good.names[1]
    }
    if (log.scale) {
      boxplot(log(y) ~ x,
              data = aux.data.to.plot, col = "lavender", xlab = x.label, main = stratum,
              ylab = paste("log(", good.names[2], ")", collapse = ""), ylim = y.limits
      )
      stripchart(log(y) ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(y ~ x,
              data = aux.data.to.plot, col = "lavender", xlab = x.label, main = stratum,
              ylab = good.names[2], ylim = y.limits
      )
      stripchart(y ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    }
    plotted <- TRUE
  }
  
  if (is.element(type.1, c("categorical", "binary")) &
      is.element(type.2, c("categorical", "binary"))) {
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    contingency.table <- table(data.to.plot)
    frequency.table <- contingency.table / rowSums(contingency.table)
    barplot(t(frequency.table),
            beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
            ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)), main = stratum,
            xlab = good.names[1]
    )
    legend(
      x = "top", legend = dimnames(frequency.table)$y, horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:ncol(frequency.table)], title = good.names[2]
    )
    plotted <- TRUE
  }
  
  if (is.element(type.1, c("categorical", "binary", "ordinal")) & type.2 == "continuous") {
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    if (log.scale) {
      boxplot(log(y) ~ x,
              data = data.to.plot, col = "lavender", xlab = good.names[1], main = stratum,
              ylab = paste("log(", good.names[2], ")", collapse = ""), ylim = y.limits
      )
      stripchart(log(y) ~ x, data = data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(y ~ x,
              data = data.to.plot, col = "lavender", xlab = good.names[1], main = stratum,
              ylab = good.names[2], ylim = y.limits
      )
    }
    plotted <- TRUE
  }
  
  if (is.element(type.2, c("categorical", "binary", "ordinal")) & type.1 == "continuous") {
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    if (log.scale) {
      boxplot(log(x) ~ y,
              data = data.to.plot, col = "lavender", xlab = good.names[2], main = stratum,
              ylab = paste("log(", good.names[1], ")", collapse = ""), ylim = x.limits
      )
      stripchart(log(x) ~ y, data = data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(x ~ y,
              data = data.to.plot, col = "lavender", xlab = good.names[2], main = stratum,
              ylab = good.names[1], ylim = x.limits
      )
      stripchart(x ~ y, data = data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    }
    plotted <- TRUE
  }
  
  if ((type.1 == "ordinal" & type.2 == "ordinal") |
      (type.1 == "binary" & type.2 == "ordinal") | (type.1 == "categorical" & type.2 == "ordinal")) {
    aux.table.1 <- table(data.to.plot[, 1])
    aux.table.2 <- table(data.to.plot[, 2])
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    if (length(aux.table.2) > length(aux.table.1)) {
      aux.means <- aggregate(y ~ x, data = data.to.plot, mean)
      boxplot(y ~ x,
              data = data.to.plot, col = "lavender", xlab = good.names[1], main = stratum,
              ylab = good.names[2], main = stratum, ylim = y.limits
      )
      points(y ~ as.factor(x),
             data = aux.means, col = "red", xlab = good.names[1], ylab = good.names[2],
             cex = 1.25, pch = 19
      )
    } else {
      aux.means <- aggregate(x ~ y, data = data.to.plot, mean)
      boxplot(x ~ y,
              data = data.to.plot, col = "lavender", xlab = good.names[2], main = stratum,
              ylab = good.names[1], main = stratum, ylim = x.limits
      )
      points(x ~ as.factor(y),
             data = aux.means, col = "red", pch = 19, cex = 1.25,
             xlab = good.names[2], ylab = good.names[1]
      )
    }
    plotted <- TRUE
  }
  
  if (type.1 == "continuous" & type.2 == "continuous") {
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    plot(y ~ x,
         data = data.to.plot, col = "blue", xlab = good.names[1], ylab = good.names[2],
         main = stratum, ylim = y.limits, xlim = x.limits
    )
    plotted <- TRUE
  }
  
  if (type.1 == "ordinal" & type.2 == "dirac.and.continuous") {
    par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    aux.table <- table(data.to.plot[, 2])
    atom <- as.numeric(dimnames(aux.table)[[1]][which.max(as.vector(aux.table))])
    aux.data.to.plot <- base::transform(data.to.plot,
                                        y = ifelse(y > atom, paste(">", atom, sep = ""),
                                                   ifelse(y < atom, paste("<", atom, sep = ""), paste("=", atom, sep = ""))
                                        )
    )
    contingency.table <- t(table(aux.data.to.plot))
    frequency.table <- contingency.table / rowSums(contingency.table)
    barplot(t(frequency.table),
            beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
            ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)),
            xlab = good.names[2]
    )
    legend(
      x = "top", legend = colnames(frequency.table), horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:ncol(frequency.table)], title = good.names[1]
    )
    aux.data.to.plot <- subset(data.to.plot, y != atom)
    if (log.scale) {
      boxplot(log(y) ~ x,
              data = aux.data.to.plot, col = "lavender", xlab = good.names[1], main = stratum,
              ylab = paste("log(", good.names[2], ")", collapse = ""), ylim = y.limits
      )
      stripchart(log(y) ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(y ~ x,
              data = aux.data.to.plot, col = "lavender", xlab = good.names[1], main = stratum,
              ylab = good.names[2], ylim = y.limits
      )
      stripchart(y ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    }
    plotted <- TRUE
  }
  
  if (type.2 == "ordinal" & type.1 == "dirac.and.continuous") {
    par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    aux.table <- table(data.to.plot[, 1])
    atom <- as.numeric(dimnames(aux.table)[[1]][which.max(as.vector(aux.table))])
    aux.data.to.plot <- base::transform(data.to.plot,
                                        x = ifelse(x > atom, paste(">", atom, sep = ""),
                                                   ifelse(x < atom, paste("<", atom, sep = ""), paste("=", atom, sep = ""))
                                        )
    )
    contingency.table <- table(aux.data.to.plot)
    frequency.table <- contingency.table / rowSums(contingency.table)
    barplot(t(frequency.table),
            beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
            ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)), main = stratum,
            xlab = good.names[1]
    )
    legend(
      x = "top", legend = colnames(frequency.table), horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:ncol(frequency.table)], title = good.names[2]
    )
    aux.data.to.plot <- subset(data.to.plot, x != atom)
    if (log.scale) {
      boxplot(log(x) ~ y,
              data = aux.data.to.plot, col = "lavender", xlab = good.names[2], main = stratum,
              ylab = paste("log(", good.names[1], ")", collapse = ""), ylim = x.limits
      )
      stripchart(log(x) ~ y, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(x ~ y,
              data = aux.data.to.plot, col = "lavender", xlab = good.names[2], main = stratum,
              ylab = good.names[1], ylim = x.limits
      )
      stripchart(x ~ y, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    }
    plotted <- TRUE
  }
  
  if (type.1 == "dirac.and.continuous" & type.2 == "dirac.and.continuous") {
    aux.table.1 <- table(data.to.plot[, 1])
    aux.table.2 <- table(data.to.plot[, 2])
    atom.1 <- as.numeric(dimnames(aux.table.1)[[1]][which.max(as.vector(aux.table.1))])
    atom.2 <- as.numeric(dimnames(aux.table.2)[[1]][which.max(as.vector(aux.table.2))])
    aux.data.to.plot <- base::transform(data.to.plot,
                                        x = ifelse(x > atom.1, paste(">", atom.1, sep = ""),
                                                   ifelse(x < atom.1, paste("<", atom.1, sep = ""), paste("=", atom.1, sep = ""))
                                        ),
                                        y = ifelse(y > atom.2, paste(">", atom.2, sep = ""),
                                                   ifelse(y < atom.2, paste("<", atom.2, sep = ""), paste("=", atom.2, sep = ""))
                                        )
    )
    contingency.table <- table(aux.data.to.plot)
    frequency.table <- contingency.table / rowSums(contingency.table)
    aux.data.to.plot <- subset(data.to.plot, x != atom.1 & y != atom.2)
    if (nrow(aux.data.to.plot) >= 1) {
      par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    }
    barplot(t(frequency.table),
            beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
            ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)), main = stratum,
            xlab = good.names[1]
    )
    legend(
      x = "top", legend = colnames(frequency.table), horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:ncol(frequency.table)], title = good.names[2]
    )
    if (nrow(aux.data.to.plot) >= 1) {
      plot(y ~ x,
           data = aux.data.to.plot, col = "blue", xlab = good.names[1],
           ylab = good.names[2], main = stratum, ylim = y.limits, xlim = x.limits
      )
      stripchart(y ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    }
    plotted <- TRUE
  }
  
  if (type.1 == "dirac.and.continuous" & type.2 == "continuous") {
    par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    aux.table <- table(data.to.plot[, 1])
    atom <- as.numeric(dimnames(aux.table)[[1]][which.max(as.vector(aux.table))])
    aux.data.to.plot <- base::transform(data.to.plot,
                                        x = ifelse(x > atom, paste(">", atom, sep = ""),
                                                   ifelse(x < atom, paste("<", atom, sep = ""), paste("=", atom, sep = ""))
                                        )
    )
    if (log.scale) {
      boxplot(log(y) ~ x,
              data = aux.data.to.plot, col = "lavender", xlab = good.names[1], main = stratum,
              ylab = paste("log(", good.names[2], ")", collapse = ""), ylim = y.limits
      )
      stripchart(log(y) ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(y ~ x,
              data = aux.data.to.plot, col = "lavender", xlab = good.names[1], main = stratum,
              ylab = good.names[2], ylim = y.limits
      )
      stripchart(y ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
      aux.data.to.plot <- subset(data.to.plot, x != atom)
      plot(y ~ x,
           data = aux.data.to.plot, col = "blue", xlab = good.names[1], ylab = good.names[2],
           main = stratum, ylim = y.limits, xlim = x.limits
      )
    }
    plotted <- TRUE
  }
  
  if (type.2 == "dirac.and.continuous" & type.1 == "continuous") {
    par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    aux.table <- table(data.to.plot[, 2])
    atom <- as.numeric(dimnames(aux.table)[[1]][which.max(as.vector(aux.table))])
    aux.data.to.plot <- base::transform(data.to.plot,
                                        y = ifelse(y > atom, paste(">", atom, sep = ""),
                                                   ifelse(y < atom, paste("<", atom, sep = ""), paste("=", atom, sep = ""))
                                        )
    )
    if (log.scale) {
      boxplot(log(x) ~ y,
              data = aux.data.to.plot, col = "lavender", xlab = good.names[2], main = stratum,
              ylab = paste("log(", good.names[1], ")", collapse = ""), ylim = x.limits
      )
      stripchart(log(x) ~ y, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(x ~ y,
              data = aux.data.to.plot, col = "lavender", xlab = good.names[2], main = stratum,
              ylab = good.names[1], ylim = y.limits
      )
      stripchart(x ~ y, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
      aux.data.to.plot <- subset(data.to.plot, y != atom)
      plot(x ~ y,
           data = aux.data.to.plot, col = "blue", xlab = good.names[2], ylab = good.names[1],
           main = stratum, xlim = x.limits, ylim = y.limits
      )
    }
    plotted <- TRUE
  }
  
  if (((type.1 == "binary" & type.2 == "ordinal") | (type.1 == "categorical" & type.2 == "ordinal")) &
      (plotted == FALSE)) {
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    contingency.table <- table(data.to.plot)
    frequency.table <- contingency.table / rowSums(contingency.table)
    barplot(t(frequency.table),
            beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
            ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)), main = stratum,
            xlab = good.names[1]
    )
    legend(
      x = "top", legend = dimnames(frequency.table)$y, horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:ncol(frequency.table)], title = good.names[2]
    )
    plotted <- TRUE
  }
  
  if (plotted == FALSE) {
    print(paste("Did not plot ", paste(columns, collapse = ","), "!", sep = ""))
  }
}



##############################################################################3

#' @title Check for taxa-metadata associations
#' @details The associations between taxa and metadata variables are tested. These are done using
#' the method/appraoched described by JA Ferreira and S Fuentes \link{https://doi.org/10.1093/bib/bbz077}
#' @param data.for.testing A data frame combined from transform_to_fractions and prep_metadata.
#' @param B Resampling default to 100000
#' @param no.of.factors Number of factor from `data.for.testing` to test.
#' @param aux.columns specific column numbers with taxa abundances.
#' @param log.scale Logical. Plot scale. Default and suggested TRUE
#' @param nominal.bound.on.FDR cut-off for multiple testing. Default=0.25
#' @param phen.data dataframe output from prep_metadata()
#' @param compare.label label to save files.
#' @param aux.TYPES How to treat each of the variables in the data.for.testing
#' @param cat.types Categorical to plot proportions positive
#' @param path_loc Location to store/save output
#' @param name.of.stratification label of stratum e.g."Sex and age-group"
#' @param verbose
#' @export
#' @examples

check_association <- function(data.for.testing,
                              B = 100000,
                              no.of.factors= NULL ,
                              aux.columns = NULL,
                              log.scale = TRUE,
                              nominal.bound.on.FDR = 0.25,
                              phen.data = NULL,
                              compare.label="case_control",
                              aux.TYPES=NULL,
                              cat.types = NULL,
                              path_loc=".",
                              name.of.stratification= "Sex and age-group",
                              verbose=TRUE){
  
  if(is.null(no.of.factors) | is.null(aux.columns) | is.null(phen.data) | is.null(aux.TYPES) | is.null(cat.types)){
    warning("Please specify no.of.factors, aux.columns, phen.data, aux.TYPES, cat.types")
    stop("These arguments cannot be NULL check on eor more of the arguments")
    
  }
  
  frequency.of.strata <- table(phen.data$stratum)
  variables.of.interest <- names(data.for.testing)
  name.of.stratification <- name.of.stratification
  for(f in (1:no.of.factors)){
    
    start.time <- Sys.time()
    
    columns.f <- cbind(f,aux.columns)
    
    #head(columns.f)
    tests <- NULL
    cat.types <- cat.types
    aux.TYPES <- aux.TYPES
    #if(f==14){
    #  aux.TYPES[42:942] <- rep("dirac.and.continuous",901)
    #}
    # STEP 1
    for(i in (1:nrow(columns.f))){
      
      test.i <- test.ONE.association.0(as.vector(columns.f[i,]),data.for.testing,
                                       aux.TYPES,variables.of.interest)
      
      if(verbose==TRUE){
        message(paste0("Processing.. ", test.i[1], " vs ", test.i[2]))
      }
      
      p.value.i <- as.numeric(test.i[3])
      #p.value.i
      if(is.na(p.value.i)){p.value.i <- 2; B.i <- 2*B}
      if(p.value.i<1 & p.value.i>0){
        B.i <- min(max(B,ceiling(400/((1-p.value.i)*p.value.i))),1000000*10*10*10)
      }else{
        B.i <- B
      }
      
      if(B.i<=B){
        test.i <- test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,B.i)
        p.value.i <- as.numeric(test.i[3])
        no.of.successes <- p.value.i*B.i
        CI.p.value.i <- Wilson.interval(no.of.successes,B.i,0.95)
        CI.p.value.i <- paste("[",CI.p.value.i[1],",",CI.p.value.i[2],"]",collapse="")
      }else{
        #print(c(i,B.i))
        tests.i <- test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))
        
        p.values.i <- as.numeric(tests.i[,3])
        no.of.successes <- sum(p.values.i*ceiling(B.i/10))
        CI.p.value.i <- Wilson.interval(no.of.successes,10*ceiling(B.i/10),0.95)
        p.value.i <- mean(p.values.i)
        CI.p.value.i <- paste("[",CI.p.value.i[1],",",CI.p.value.i[2],"]",collapse="")
        
      }
      
      output.i <- data.frame(i=as.vector(columns.f[i,])[1],j=as.vector(columns.f[i,])[2],
                             variable.1=as.character(test.i[1]),variable.2=as.character(test.i[2]),
                             p.value=p.value.i,
                             sign=as.numeric(test.i[4]),
                             test=as.character(test.i[5]),
                             sampling.characteristics=as.character(test.i[6]),
                             CI.p.value=CI.p.value.i)
      #print(output.i)
      tests <- rbind(tests,output.i)
      #tests
      
    }
    
    aux.tests <- tests; rownames(aux.tests) <- NULL
    aux.tests <- as.data.frame(aux.tests)
    names(aux.tests) <- c("i","j","variable.1","variable.2","p.value","sign","test","sampling.characteristics","CI.95.pc.p.value")
    
    aux.tests$i <- as.numeric(as.character(aux.tests$i))
    aux.tests$j <- as.numeric(as.character(aux.tests$j))
    aux.tests$variable.1 <- as.character(aux.tests$variable.1)
    aux.tests$variable.2 <- as.character(aux.tests$variable.2)
    aux.tests$p.value <- as.numeric(as.character(aux.tests$p.value))
    aux.tests$sign <- as.numeric(as.character(aux.tests$sign))
    aux.tests$test <- as.character(aux.tests$test)
    aux.tests$sampling.characteristics <- as.character(aux.tests$sampling.characteristics)
    #head(aux.tests); str(aux.tests)
    
    path.save.assoc <- paste0(path_loc, "associations.between.",compare.label,".",names(data.for.testing)[f],".and.bacteria.CSV")
    
    #file.name <- paste("ms/taxa_test/associations.between.", compare.label,".", names(data.for.testing)[f],"and.bacteria.CSV",sep=".")
    write.csv2(aux.tests,file=path.save.assoc,row.names=FALSE,quote=TRUE)
    
    multiple.tests <- na.omit(aux.tests[,-ncol(aux.tests)])
    multiple.tests$variable.1 <- as.character(multiple.tests$variable.1)
    multiple.tests$variable.2 <- as.character(multiple.tests$variable.2)
    multiple.tests$test <- as.character(multiple.tests$test)
    multiple.tests <- multiple.tests[order(multiple.tests$p.value),]; row.names(multiple.tests) <- NULL
    aux.histogram <- hist(multiple.tests$p.value,plot=FALSE)
    gamma.hat <- min(1,aux.histogram$density[length(aux.histogram$density)]); gamma.hat
    multiple.tests <- base::transform(multiple.tests,rank=(1:nrow(multiple.tests)))
    bound.FDR <- formatC(nrow(multiple.tests)*multiple.tests$p.value/(1:nrow(multiple.tests)),
                         digits=4,format="f")
    better.bound.FDR <- formatC(gamma.hat*nrow(multiple.tests)*multiple.tests$p.value/(1:nrow(multiple.tests)),
                                digits=4,format="f")
    multiple.tests <- base::transform(multiple.tests,
                                      bound.FDR=as.numeric(bound.FDR),
                                      better.bound.FDR=as.numeric(better.bound.FDR))
    #head(multiple.tests); str(multiple.tests)
    
    if(any(multiple.tests$better.bound.FDR<=nominal.bound.on.FDR)){
      
      cut.off <- max((1:nrow(multiple.tests))[multiple.tests$better.bound.FDR<=nominal.bound.on.FDR])
      selected.tests <- multiple.tests[1:cut.off,]
      rownames(selected.tests) <- NULL
      path.save.select <- paste0(path_loc, "SELECTED.associations.between.",compare.label,".",names(data.for.testing)[f],".and.bacteria.CSV")
      
      #path.save <- paste0("ms/taxa_test/iliv1v2/", "SELECTED.associations.between.",compare.label,".",names(data.for.testing)[f],".and.bacteria.CSV",sep=".")
      #file.name <- paste("ms/taxa_test/SELECTED.associations.between", compare.label,".",names(data.for.testing)[f],"and.bacteria.CSV",sep=".")
      write.csv2(selected.tests,file=path.save.select,row.names=FALSE,quote=TRUE)
      path.save.select2 <- paste0(path_loc, "SELECTED.associations.between.",compare.label,".",names(data.for.testing)[f],".and.bacteria.pdf")
      #file.name <- paste("ms/taxa_test/SELECTED.associations.between", compare.label,".",names(data.for.testing)[f],"and.bacteria.pdf",sep=".")
      pdf(file=path.save.select2,width=9,height=7)
      
      par(mfrow=c(1,1))
      hist(multiple.tests$p.value,col="lavender",xlab="P-value",
           main=paste("Associations between ",names(data.for.testing)[f],
                      " and bacterial abundance;\n estimated proportion of potential associations: ",
                      formatC(1-gamma.hat,digits=2,format="f"),sep=""),
           prob=TRUE,cex.lab=1,cex.axis=1,cex.main=1)
      abline(h=1,col="red")
      abline(h=gamma.hat,col="blue",lty=2)
      
      for(r in (1:nrow(selected.tests))){
        columns.r <- as.numeric(selected.tests[r,1:2])
        illustrate.ONE.association(columns=columns.r,data.for.testing,cat.types,log.scale)
      }
      
      aux.data.for.testing <- subset(data.for.testing,
                                     select=c(selected.tests[1,3],selected.tests[,4],"stratum"))
      aux.types <- cat.types[c(selected.tests[1,1],selected.tests[,2])]
      frequency.of.strata <- table(aux.data.for.testing$stratum)
      if(is.factor(aux.data.for.testing[,1]) | is.character(aux.data.for.testing[,1])){
        min.by.stratum <- aggregate(as.formula(paste(names(aux.data.for.testing)[1],"stratum",sep="~")),
                                    data=aux.data.for.testing,min)
        max.by.stratum <- aggregate(as.formula(paste(names(aux.data.for.testing)[1],"stratum",sep="~")),
                                    data=aux.data.for.testing,max)
        good.strata <- min.by.stratum[max.by.stratum[,2]!=min.by.stratum[,2],1]
      }else{
        sd.by.stratum <- aggregate(as.formula(paste(names(aux.data.for.testing)[1],"stratum",sep="~")),
                                   data=aux.data.for.testing,sd)
        good.strata <- sd.by.stratum[sd.by.stratum[,2]>0,1]
      }
      
      frequency.of.strata <- frequency.of.strata[is.element(dimnames(frequency.of.strata)[[1]],good.strata)]
      
      minimal.sample.size <- 25
      strata.of.interest <- dimnames(frequency.of.strata)[[1]][frequency.of.strata>=minimal.sample.size]
      strata.of.interest
      
      for(r in (2:(ncol(aux.data.for.testing)-1))){
        
        if(log.scale){
          x.limits <- range(aux.data.for.testing[,1],na.rm=TRUE)
          y.limits <- range(log(aux.data.for.testing[aux.data.for.testing[,r]>0,r]),na.rm=TRUE)
        }else{
          x.limits <- range(aux.data.for.testing[,1],na.rm=TRUE)
          y.limits <- range(aux.data.for.testing[,r],na.rm=TRUE)
        }
        
        for(s in strata.of.interest){
          
          data.for.testing.s <- na.omit(subset(aux.data.for.testing,stratum==s)[c(1,r)])
          title.s <- paste(name.of.stratification," stratum: ",s,sep="")
          if(nrow(data.for.testing.s)>=minimal.sample.size & sd(data.for.testing.s[,2])>0){
            illustrate.ONE.association.BY.STRATUM(data.for.testing.s,aux.types[c(1,r)],stratum=title.s,
                                                  x.limits,y.limits,log.scale=TRUE)
          }
        }
      }
      
      dev.off()
      
    }
    
    end.time <- Sys.time()
    time.taken <- end.time-start.time
    print(time.taken)
    
  }
  
}


##############################################################################3
path.files <- "ms/taxa_test/ctrl_ili/"
pattern.name <- "SELECTED.associations.between.ctrl_ili."
save_xlsx_workbook <- function(path.files=NULL, 
                               pattern.name="SELECTED.associations.between.ctrl_ili."){
  require(xlsx)
  files_list <- list.files(path.files, pattern = pattern.name)
  dbf.files <- files_list[-grep(".pdf", files_list, fixed=T)]
  dbf.files.path <- paste0(path.files,dbf.files)
  
  # create first excel workbook
  file1 <- read_delim(dbf.files.path[1], ";", 
                      escape_double = FALSE, 
                      trim_ws = TRUE)
  
  y <- dbf.files.path[1]
  y <- gsub(pattern.name, "", y)
  y <- gsub(path.files, "", y)
  y <- gsub(".and.bacteria.CSV", "", y)
  
  file.name.to.save <- paste0(path.files,"Taxa_Associations.xlsx")
  write.xlsx(file1, file=file.name.to.save,
             sheetName=y, append=FALSE)
  
  for(i in dbf.files.path[-1]){
    
    x <- dbf.files[i]
    x <- gsub(pattern.name,
              "", i)
    x <- gsub(path.files, "", x)
    x <- gsub(".and.bacteria.CSV", "", x)
    #print(x)
    
    filex <- read_delim(i, ";", 
                        escape_double = FALSE, 
                        trim_ws = TRUE)
    
    write.xlsx(filex, 
               file=file.name.to.save,
               sheetName=x, append=TRUE) 
  }
  
}


##############################################################################3

combine_tables <- function(path.files,pattern.name){
  files_list <- dbf.files <- dbf.files.path <- all_test <- all_test_df <- NULL
  files_list <- list.files(path.files, pattern = pattern.name)
  dbf.files <- files_list[-grep(".pdf", files_list, fixed=T)]
  dbf.files.path <- paste0(path.files,dbf.files)
  all_test <- NULL
  for (j in dbf.files.path){
    # create first excel workbook
    file1 <- read_delim(j, ";", 
                        escape_double = FALSE, 
                        trim_ws = TRUE)
    
    y <- j
    y <- gsub(pattern.name, "", y)
    y <- gsub(path.files, "", y)
    y <- gsub(".and.bacteria.CSV", "", y)
    file1$variable <- y
    all_test <- rbind(all_test,file1)
  }
  
  all_test_df <- as.data.frame(all_test)
  all_test_df$p.value <- as.numeric(gsub(",", ".", all_test_df$p.value))
  all_test_df$bound.FDR <- as.numeric(gsub(",", ".", all_test_df$bound.FDR))
  all_test_df$better.bound.FDR <- as.numeric(gsub(",", ".", all_test_df$better.bound.FDR))
  
  colnames(all_test_df) <- c("var.i","var.j", "variable", "taxa", 
                             "p.value", "sign", "test", "sampling.characteristics",
                             "rank.within.variable.tested", "bound.FDR", 
                             "better.bound.FDR","variable2")
  return(all_test_df)
}

###############################################################################
plot_taxa_associations <- function(tab, ps, 
                                   adj.pval=TRUE, 
                                   group="condition_status",
                                   cut.off=0.25,
                                   p.widths = c(2,1)) {
  
  require(patchwork)
  require(tibble)
  tab$variable <- gsub("_2014","", tab$variable)
  tab$taxa <- gsub("Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium", 
                   "A-N-P-R*",	tab$taxa)
  
  #group <- "condition_status"
  # calculate prev
  ps.ctrl <- prevalence(subset_samples(ps, condition_status=="control"))
  ps.ctrl <- as.data.frame(ps.ctrl) %>% 
    rownames_to_column("taxa") %>% mutate(group="control")
  colnames(ps.ctrl) <- c("taxa", "prevalence", "group")
  ps.ili <- prevalence(subset_samples(ps, condition_status=="acute_ili"))
  ps.ili <- as.data.frame(ps.ili) %>% 
    rownames_to_column("taxa") %>% mutate(group="acute_ili")
  colnames(ps.ili) <- c("taxa", "prevalence", "group")
  
  ps.prev<- bind_rows(ps.ili, ps.ctrl)
  ps.prev$taxa <- gsub("Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium", 
                       "A-N-P-R*",	ps.prev$taxa)
  
  ps.prev
  
  #ps.prev <- ps.prev %>% 
  #filter(taxa %in% unique(tab$taxa)) %>% 
  # arrange(desc(prevalence))
  
  
  ##############################plot associations########################################
  # plot associations 
  tab2 <- tab %>% 
    filter(better.bound.FDR <= cut.off)
  
  ps.prev <- ps.prev %>% 
    filter(taxa %in% unique(tab2$taxa)) 
  
  order_tx <- ps.prev %>% 
    filter(group=="acute_ili") %>% 
    arrange(desc(prevalence))
  
  tab2$taxa <- factor(tab2$taxa, levels = unique(order_tx$taxa))
  tab2$Association <- ifelse(tab2$sign == -1, "Negative", "Positive")
  
  p.associations <- ggplot(tab2, aes(variable, taxa)) +
    geom_tile(aes(fill= Association),size=0.15) +
    #geom_text(aes(label = sign), size=2, color="white") +
    scale_fill_manual(values = c(Positive = "#457b9d", Negative = "#e63946" )) +
    theme_minimal(base_size = 8) + 
    theme(axis.text.x = element_text(angle=90, 
                                     hjust=0.95,
                                     vjust=0.2),
          axis.text.y = element_text(face = "italic"),
          axis.title = element_blank(),
          legend.position = "top",
          plot.caption = element_text(face = "italic", hjust = -1.0)) 
  #+    labs(caption = "*A-N-P-R : Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium")
  
  ##############################plot prev########################################
  # create a gradient color palette for prevalence
  gray_grad <- colorRampPalette(c("white","#161a1d"))
  gray_grad_cols <- gray_grad(5)
  
  ps.prev$taxa <- factor(ps.prev$taxa, levels = unique(order_tx$taxa))
  
  p.prv <- ggplot(ps.prev, aes(group, taxa))+ 
    #geom_point(aes(fill=mean_abundance, size=sd_abundance), shape=21) +
    geom_tile(aes(fill= prevalence),colour="grey50",size=0.25) +
    #geom_text(aes(label = round(sd_abundance, 3)), size=2,color="white") +
    scale_fill_gradientn(colours = gray_grad_cols)+
    theme_void(base_size = 8) + 
    theme(axis.text.x = element_text(angle=90, 
                                     hjust=0.95,
                                     vjust=0.2),
          axis.title = element_blank(),
          legend.position = "right",
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(fill=NA))
  
  p.comb <- p.associations + p.prv + plot_layout(widths = p.widths)
  return(p.comb) 
  
}


##########################################3
plot_associations2 <- function(tab, ps, adj.pval=TRUE, 
                               group="condition_status",
                               cut.off=0.25,
                               p.widths = c(2,1)) {
  
  require(patchwork)
  
  tab$variable <- gsub("_2014","", tab$variable)
  tab$taxa <- gsub("Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium", 
                   "A-N-P-R*",	tab$taxa)
  
  #group <- "condition_status"
  # calculate prev
  ps.ctrl <- prevalence(subset_samples(ps, condition_status=="acute_ili_14_days"))
  ps.ctrl <- as.data.frame(ps.ctrl) %>% 
    rownames_to_column("taxa") %>% mutate(group="acute_ili_14_days")
  colnames(ps.ctrl) <- c("taxa", "prevalence", "group")
  ps.ili <- prevalence(subset_samples(ps, condition_status=="acute_ili"))
  ps.ili <- as.data.frame(ps.ili) %>% 
    rownames_to_column("taxa") %>% mutate(group="acute_ili")
  colnames(ps.ili) <- c("taxa", "prevalence", "group")
  
  ps.prev<- bind_rows(ps.ili, ps.ctrl)
  ps.prev$taxa <- gsub("Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium", 
                       "A-N-P-R*",	ps.prev$taxa)
  #ps.prev <- ps.prev %>% 
  #filter(taxa %in% unique(tab$taxa)) %>% 
  # arrange(desc(prevalence))
  
  
  ##############################plot associations########################################
  # plot associations 
  if(adj.pval==TRUE) {
    
    tab2 <- tab %>% 
      filter(better.bound.FDR <= cut.off)
    
    #ps.prev <- ps.prev %>% 
    #  filter(taxa %in% unique(tab2$taxa)) %>% 
    #  arrange(desc(prevalence))
    
    #tab2$taxa <- factor(tab2$taxa, levels = unique(ps.prev$taxa))
    p.associations <- ggplot(tab2, aes(variable, taxa)) 
    p.associations <- p.associations + 
      geom_tile(aes(fill= log10(better.bound.FDR)),colour="white",size=0.15) 
  } else{
    
    tab2 <- tab %>% 
      filter(p.value <= cut.off)
    
    #ps.prev <- ps.prev %>% 
    #   filter(taxa %in% unique(tab2$taxa))%>% 
    # arrange(desc(prevalence))
    
    #tab2$taxa <- factor(tab2$taxa, levels = unique(ps.prev$taxa))
    
    p.associations <- ggplot(tab2, aes(variable, taxa)) 
    p.associations <- p.associations + 
      geom_tile(aes(fill= log10(p.value)),colour="white",size=0.15) 
  }
  
  p.associations <- p.associations + 
    geom_text(aes(label = sign), size=2, color="white") +
    scale_fill_gradient(low = "#e63946", high ="#457b9d")+
    theme_biome_utils() + 
    theme(axis.text.x = element_text(angle=90, 
                                     hjust=0.95,
                                     vjust=0.2),
          axis.text.y = element_text(face = "italic"),
          axis.title = element_blank(),
          legend.position = "top",
          plot.caption = element_text(face = "italic", hjust = -1.0)) 
  #+    labs(caption = "*A-N-P-R : Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium")
  
  ##############################plot prev########################################
  # create a gradient color palette for prevalence
  gray_grad <- colorRampPalette(c("white","#161a1d"))
  gray_grad_cols <- gray_grad(5)
  ps.prev <- ps.prev %>% 
    filter(taxa %in% unique(tab2$taxa))
  p.prv <- ggplot(ps.prev, aes(group, taxa))+ 
    #geom_point(aes(fill=mean_abundance, size=sd_abundance), shape=21) +
    geom_tile(aes(fill= prevalence),colour="white",size=0.15) +
    #geom_text(aes(label = round(sd_abundance, 3)), size=2,color="white") +
    scale_fill_gradientn(colours = gray_grad_cols)+
    theme_void() + 
    theme(axis.text.x = element_text(angle=90, 
                                     hjust=0.95,
                                     vjust=0.2),
          axis.title = element_blank(),
          legend.position = "right",
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank())
  
  p.comb <- p.associations + p.prv + plot_layout(widths = p.widths)
  return(p.comb) 
  
}

#################################################################################

plot_associations3 <- function(tab, ps, adj.pval=TRUE, 
                               group="condition_status",
                               cut.off=0.25,
                               p.widths = c(2,1)) {
  
  require(patchwork)
  
  tab$variable <- gsub("_2014","", tab$variable)
  tab$taxa <- gsub("Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium", 
                   "A-N-P-R*",	tab$taxa)
  
  #group <- "condition_status"
  # calculate prev
  ps.ctrl <- prevalence(subset_samples(ps, condition_status=="acute_ili_recovery"))
  ps.ctrl <- as.data.frame(ps.ctrl) %>% 
    rownames_to_column("taxa") %>% mutate(group="acute_ili_recovery")
  colnames(ps.ctrl) <- c("taxa", "prevalence", "group")
  ps.ili <- prevalence(subset_samples(ps, condition_status=="acute_ili"))
  ps.ili <- as.data.frame(ps.ili) %>% 
    rownames_to_column("taxa") %>% mutate(group="acute_ili")
  colnames(ps.ili) <- c("taxa", "prevalence", "group")
  
  ps.prev<- bind_rows(ps.ili, ps.ctrl)
  ps.prev$taxa <- gsub("Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium", 
                       "A-N-P-R*",	ps.prev$taxa)
  #ps.prev <- ps.prev %>% 
  #filter(taxa %in% unique(tab$taxa)) %>% 
  # arrange(desc(prevalence))
  
  
  ##############################plot associations########################################
  # plot associations 
  if(adj.pval==TRUE) {
    
    tab2 <- tab %>% 
      filter(better.bound.FDR <= cut.off)
    
    #ps.prev <- ps.prev %>% 
    #  filter(taxa %in% unique(tab2$taxa)) %>% 
    #  arrange(desc(prevalence))
    
    #tab2$taxa <- factor(tab2$taxa, levels = unique(ps.prev$taxa))
    p.associations <- ggplot(tab2, aes(variable, taxa)) 
    p.associations <- p.associations + 
      geom_tile(aes(fill= log10(better.bound.FDR)),colour="white",size=0.15) 
  } else{
    
    tab2 <- tab %>% 
      filter(p.value <= cut.off)
    
    #ps.prev <- ps.prev %>% 
    #   filter(taxa %in% unique(tab2$taxa))%>% 
    # arrange(desc(prevalence))
    
    #tab2$taxa <- factor(tab2$taxa, levels = unique(ps.prev$taxa))
    
    p.associations <- ggplot(tab2, aes(variable, taxa)) 
    p.associations <- p.associations + 
      geom_tile(aes(fill= log10(p.value)),colour="white",size=0.15) 
  }
  
  p.associations <- p.associations + 
    geom_text(aes(label = sign), size=2, color="white") +
    scale_fill_gradient(low = "#e63946", high ="#457b9d")+
    theme_biome_utils() + 
    theme(axis.text.x = element_text(angle=90, 
                                     hjust=0.95,
                                     vjust=0.2),
          axis.text.y = element_text(face = "italic"),
          axis.title = element_blank(),
          legend.position = "top",
          plot.caption = element_text(face = "italic", hjust = -1.0)) 
  #+    labs(caption = "*A-N-P-R : Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium")
  
  ##############################plot prev########################################
  # create a gradient color palette for prevalence
  gray_grad <- colorRampPalette(c("white","#161a1d"))
  gray_grad_cols <- gray_grad(5)
  ps.prev <- ps.prev %>% 
    filter(taxa %in% unique(tab2$taxa))
  p.prv <- ggplot(ps.prev, aes(group, taxa))+ 
    #geom_point(aes(fill=mean_abundance, size=sd_abundance), shape=21) +
    geom_tile(aes(fill= prevalence),colour="white",size=0.15) +
    #geom_text(aes(label = round(sd_abundance, 3)), size=2,color="white") +
    scale_fill_gradientn(colours = gray_grad_cols)+
    theme_void() + 
    theme(axis.text.x = element_text(angle=90, 
                                     hjust=0.95,
                                     vjust=0.2),
          axis.title = element_blank(),
          legend.position = "right",
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank())
  
  p.comb <- p.associations + p.prv + plot_layout(widths = p.widths)
  return(p.comb) 
  
}

#############################################
plot_paired <- function (x, select.taxa = NULL, group = NULL, group.colors = NULL, 
                         dot.opacity = 0.25, dot.size = 2, add.box = FALSE, box.opacity = 0.25, 
                         group.order = NULL, add.violin = TRUE, violin.opacity = 0.25, 
                         ncol = NULL, nrow = NULL, line = NULL, line.down = "#7209b7", 
                         line.stable = "#8d99ae", line.up = "#14213d", line.na.value = "grey50", 
                         line.guide = "legend", line.opacity = 0.25, line.size = 1, 
                         jitter.width = 0) 
{
  Abundance <- change <- change.order <- change.sign <- linevar <- NULL
  if (is.na(line) | is.null(line)) {
    stop(" 'line' argument cannot be empty")
  }
  if (is.na(group) | is.null(group)) {
    stop(" 'group' argument cannot be empty")
  }
  xmeta <- meta(x)
  for (i in select.taxa) {
    df.tx <- as.data.frame(abundances(x)[i, ])
    colnames(df.tx) <- i
    xmeta <- cbind(xmeta, df.tx)
  }
  #if (length(unique(xmeta[, group])) > 2) {
  #  stop("Only two group comparison e.g. before n after")
  #}
  if (!is.factor(xmeta[, group])) {
    xmeta$group <- factor(as.character(xmeta[, group]))
  }
  xmeta_lf <- xmeta %>% pivot_longer(cols = all_of(select.taxa), 
                                     names_to = "taxa", values_to = "Abundance")
  xmeta_lf$linevar <- factor(xmeta_lf[[line]])
  x.grp <- sym(group)
  df2 <- suppressWarnings(xmeta_lf %>% arrange(taxa, linevar, 
                                               !!x.grp) %>% group_by(taxa, linevar) %>% 
                            summarise(change = diff(Abundance)))
  xmeta_lf_2 <- suppressWarnings(xmeta_lf %>% arrange(taxa, 
                                                      linevar, !!x.grp) %>% group_by(taxa, linevar))
  xmeta_lf_2 <- xmeta_lf_2 %>% left_join(df2)
  xmeta_lf_2$change.sign <- sign(xmeta_lf_2$change)
  if (!is.null(group.order)) {
    xmeta_lf_2[, group] <- factor(xmeta_lf_2[, group], levels = group.order)
  }
  xmeta_lf_2 <- xmeta_lf_2 %>% mutate(change.order = ifelse(change.sign == 
                                                              1, "Up", ifelse(change.sign == -1, "Down", "Stable")))
  p <- ggplot(data = xmeta_lf_2, aes_string(x = "group", y = "Abundance", 
                                            fill = "group")) + geom_point(aes_string(x = "group", 
                                                                                     fill = "group"), position = position_jitter(width = jitter.width), 
                                                                          size = dot.size, alpha = dot.opacity, shape = 21) + geom_line(aes(group = linevar, 
                                                                                                                                            color = change.order), size = line.size, alpha = line.opacity)
  p <- p + scale_fill_manual(values = group.colors)
  if (add.box == TRUE) {
    p <- p + geom_boxplot(width = 0.2, outlier.shape = NA, 
                          alpha = box.opacity)
  }
  if (add.violin == TRUE) {
    p <- p + geom_half_violin(data = xmeta_lf_2 %>% filter(group == 
                                                             unique(xmeta_lf_2$group)[1]), position = position_nudge(x = -0.15, 
                                                                                                                     y = 0), alpha = violin.opacity, side = "l") + geom_half_violin(data = xmeta_lf_2 %>% 
                                                                                                                                                                                      filter(group == unique(xmeta_lf_2$group)[2]), position = position_nudge(x = 0.15, 
                                                                                                                                                                                                                                                              y = 0), alpha = violin.opacity, side = "r")
  }
  if (length(select.taxa) >= 2) {
    p <- p + facet_wrap(~taxa, scales = "free", ncol = ncol, 
                        nrow = nrow)
  }
  p <- p + scale_color_manual(values = c(Up = line.up, Down = line.down, 
                                         Stable = line.stable), guide = line.guide)
  p <- p + theme_biome_utils() + xlab(group)
  return(p)
}



