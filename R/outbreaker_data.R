#' Process input data for outbreaker
#'
#' This function performs various checks on input data given to outbreaker.  It
#' takes a list of named items as input, performs various checks, set defaults
#' where arguments are missing, and return a correct list of data input. If no
#' input is given, it returns the default settings.
#'
#' Acceptable arguments for ... are:
#' \describe{
#' \item{dates}{a vector indicating the collection dates, provided either as
#' integer numbers or in a usual date format such as \code{Date} or
#' \code{POSIXct} format. By convention, zero will indicate the oldest date. 
#' Cases must be ordering by ascending onset date.}
#'
#' \item{age_group}{a vector indicating the age group of the cases, provided as
#' integer numbers. The value of age group corresponds to the position of this age
#' group in a_dens.}
#' 
#' \item{region}{a vector indicating the region of the cases, provided as
#' integer numbers or characters. If numeric, the value of the region corresponds 
#' to the position of the region in the distance matrix and the population vector.
#' Otherwise, the value corresponds to the region and will be matched to the distance
#' matrix and the population vector.}
#' 
#' \item{w_dens}{a vector of numeric values indicating the generation time
#' distribution, reflecting the infectious potential of a case t = 1, 2, ...
#' time steps after infection. By convention, it is assumed that
#' newly infected patients cannot see new infections on the same time step. If not
#' standardized, this distribution is rescaled to sum to 1.}
#'
#' \item{f_dens}{similar to \code{w_dens}, except that this is the distribution
#' of the incubation period, i_e. time interval between the reported onset date
#' and the infection date.}
#' 
#' \item{a_dens}{a matrix of numeric values indicating the contact between age 
#' groups, reflecting on the infectious potential of a case for a given age group.}
#' 
#' \item{genotype}{a character vector showing the genotype in each case.}
#' 
#' \item{is_cluster}{an integer vector indicating which group of cases each case
#'  belongs to.}
#'  
#' \item{s_dens}{a matrix of numeric values indicating the initial value of 
#' the connectivity between region. Only needed if a and b are fixed in the model, 
#' otherwise NULL.}
#' 
#' \item{population}{a double vector indicating the population in every region 
#' considered in the study.}
#' 
#' \item{distance}{a double matrix indicating the distance between each region.}
#' 
#' \item{import}{a logical vector indicating whether each case is an import (TRUE)
#'  or not (FALSE).}
#'
#'}
#'
#' @param ... a list of data items to be processed (see description)
#'
#' @param data optionally, an existing list of data item as returned by \code{outbreaker_data}.
#'
#' @author Initial version by Thibaut Jombart, rewritten by Alexis Robert (\email{alexis.robert@lshtm.ac.uk})
#'
#' @import data.table
#' @export
#'
#' @return 
#' 
#' A named list containing the value of each elements listed in the 'Details' section. This list describes the data that will be used by the function \code{outbreaker()}.
#'
#' @examples
#'
#' data("toy_outbreak_short")
#' dt_cases <- toy_outbreak_short$cases
#' dt_cases <- dt_cases[order(dt_cases$Date), ]
#' dt_regions <- toy_outbreak_short$dt_regions
#' all_dist <- geosphere::distGeo(matrix(c(rep(dt_regions$long, nrow(dt_regions)), 
#'                                         rep(dt_regions$lat, nrow(dt_regions))), 
#'                                       ncol = 2), 
#'                                matrix(c(rep(dt_regions$long, each = nrow(dt_regions)), 
#'                                         rep(dt_regions$lat, each = nrow(dt_regions))),
#'                                       ncol = 2))
#' 
#' dist_mat <- matrix(all_dist/1000, nrow = nrow(dt_regions))
#' pop_vect <- dt_regions$population
#' names(pop_vect) <- rownames(dist_mat) <- colnames(dist_mat) <- dt_regions$region
#' data <- outbreaker_data(dates = dt_cases$Date, age_group = dt_cases$age_group,
#'                         region = dt_cases$Cens_tract, population = pop_vect, 
#'                         distance = dist_mat)
#'
outbreaker_data <- function(..., data = list(...)) {
  
  ## SET DEFAULTS ##
  defaults <- list(dates = NULL, region = NULL, age_group = NULL,
                   w_dens = NULL, f_dens = NULL, a_dens = NULL, N = 0L,
                   max_range = NA, can_be_ances_reg = NULL, log_w_dens = NULL, 
                   log_f_dens = NULL, log_a_dens = NULL,genotype = NULL, 
                   is_cluster = NULL, cluster = NULL, population = NULL, 
                   distance = NULL, import = NULL, s_dens = NULL, log_s_dens = NULL
                   )
  
  ## MODIFY DATA WITH ARGUMENTS ##
  data <- modify_defaults(defaults, data, FALSE)
  
  
  ## CHECK DATA ##
  ## CHECK CLUSTER
  if(is.null(data$is_cluster)){
    data$is_cluster <- rep(1, length(data$dates))
    cluster_list <- list(seq_along(data$is_cluster))
    names(cluster_list) <- 1
    data$cluster <- cluster_list
  } else{
    if(any(is.na(data$is_cluster)))
      stop("non-finite values detected in is_cluster")
    data$is_cluster <- factor(data$is_cluster, levels = unique(data$is_cluster))
    data$is_cluster <- as.numeric(data$is_cluster)
    if(all(data$is_cluster == 1)){
      cluster_list <- list(seq_along(data$is_cluster))
    } else{
      cluster_list <- sapply(unique(data$is_cluster), function(X) 
        return(which(data$is_cluster == as.numeric(X))))
    }
    names(cluster_list) <- seq_along(unique(data$is_cluster))
    data$cluster <- cluster_list
  }
  
  ## CHECK DATES
  if (!is.null(data$dates)) {
    if (inherits(data$dates, "Date")) {
      data$dates <- data$dates-min(data$dates)
    }
    if (inherits(data$dates, "POSIXct")) {
      data$dates <- difftime(data$dates, min(data$dates), units="days")
    }
    if(any(sort(data$dates) != data$dates))
      stop("Dates are not sorted, sort cases by onset dates")
    data$dates <- as.integer(round(data$dates))
    data$N <- length(data$dates)
    data$max_range <- max(vapply(unique(data$is_cluster), function(X){
      return(diff(range(data$dates[which(data$is_cluster == X)])))
    }, 1))
    ## get temporal ordering constraint:
  }
  
  ## CHECK AGE GROUP
  if (!is.null(data$age_group)) {
    if (!is.numeric(data$age_group)) {
      stop("Age group is not numeric")
    }
  }
  
  ## CHECK W_DENS
  if (!is.null(data$w_dens)) {
    if (any(data$w_dens<0)) {
      stop("w_dens has negative entries (these should be probabilities!)")
    }
    
    if (any(!is.finite(data$w_dens))) {
      stop("non-finite values detected in w_dens")
    }
    
    
    ## Remove trailing zeroes to prevent starting with -Inf temporal loglike
    if(data$w_dens[length(data$w_dens)] < 1e-15) {
      final_index <- max(which(!is.infinite(log(data$w_dens)))) - 1
      data$w_dens <- data$w_dens[1:final_index]
      warning("Removed trailing zeroes found in w_dens")
    }
    
    ## add an exponential tail to 'w'
    ## to cover the span of the outbreak
    ## (avoids starting with -Inf temporal loglike)
    if (length(data$w_dens) < data$max_range) {
      final_index <- length(data$w_dens)
      length_to_add <- (data$max_range - final_index) + 1 
      val_to_add <- min(1e-100, data$w_dens[final_index])
      data$w_dens <- c(data$w_dens[-final_index], 
                       rep(val_to_add, length_to_add))
      data$w_dens <- data$w_dens/sum(data$w_dens)
    }
    
    ## standardize the mass function
    data$w_dens <- data$w_dens / sum(data$w_dens)
    data$log_w_dens <- matrix(log(data$w_dens), nrow = 1)
  }
  
  ## CHECK F_DENS
  if (!is.null(data$w_dens) && is.null(data$f_dens)) {
    data$f_dens <- data$w_dens
  }
  if (!is.null(data$f_dens)) {
    if (any(data$f_dens<0)) {
      stop("f_dens has negative entries (these should be probabilities!)")
    }
    
    if (any(!is.finite(data$f_dens))) {
      stop("non-finite values detected in f_dens")
    }
    
    data$f_dens <- data$f_dens / sum(data$f_dens)
    data$log_f_dens <- log(data$f_dens)
  }
  
  ## CHECK A_DENS
  if (!is.null(data$a_dens)) {
    if(any(dim(data$a_dens)<max(data$age_group)))
      stop("The dimension of the contact probability matrix is lower than the maximum value of age group")
    if (any(data$a_dens<0)) {
      stop("a_dens has negative entries (these should be probabilities!)")
    }
    
    if (any(!is.finite(data$a_dens))) {
      stop("non-finite values detected in a_dens")
    }
    
    
    ## Remove trailing zeroes to prevent starting with -Inf temporal loglike
    if(any(data$a_dens < 1e-15)) {
      data$a_dens[data$a_dens<1e-15] <- 1e-15
      data$a_dens <- t(t(data$a_dens)/colSums(data$a_dens))
      warning("Removed trailing zeroes found in a_dens")
    }

    ## standardize the mass function
    data$a_dens <- t(t(data$a_dens)/colSums(data$a_dens))
    data$log_a_dens <- matrix(log(data$a_dens), nrow = nrow(data$a_dens))
    data$log_a_dens <- list(data$log_a_dens)
  }
  
  ## CHECK GENOTYPE
  if(!is.null(data$genotype)){
    data$genotype[is.na(data$genotype)] <- "Not attributed"
  } else
    data$genotype <- rep("Not attributed", 
                         length(data$dates))
  
  ## CHECK IMPORT
  if (!is.null(data$import)) {
    if(any(!is.logical(data$import)))
      stop("import should be logical")
    if(length(data$import) == 1)
      data$import <- rep(data$import, data$N) else if(length(data$import) != data$N)
        stop("import should be either of length 1 or N")
  } else {
    data$import <- rep(FALSE, data$N)
  }
  
  ## CHECK region
  if (!is.null(data$region)) {
    ## CHECK POPULATION
    if(!is.null(data$population)){
      if(any(!is.numeric(data$population)) || any(data$population<0))
        stop("The vector population can only include positive numeric values")
    } else
      stop("The population vector is null")
    
    ## CHECK DISTANCE
    if(!is.null(data$distance)){
      if(any(!is.numeric(data$distance)))
        stop("The matrix distance can only include numeric values")
      if(!is.matrix(data$distance))
        stop("distance has to be a matrix")
      if(dim(data$distance)[1] != dim(data$distance)[2])
        stop("distance has to be a square matrix")
      if(any(rownames(data$distance) != colnames(data$distance)))
        stop("The rownames and colnames of the matrix distance should be the same")
    }else
      stop("The distance matrix is null")
    if(!any(is.element(names(data$population), colnames(data$distance))))
      stop("The vector population should have the same names as the distance matrix")
    if(length(match.arg(names(data$population), 
                        colnames(data$distance), 
                        several.ok = TRUE)) != length(data$population))
      stop("The vector population should have the same names as the distance matrix")
    if(any(names(data$population) != colnames(data$distance)))
      data$distance <- data$distance[names(data$population), names(data$population)]
    
    if(!is.double(data$region)){
      if(!is.character(data$region))
        stop("region is not an integer, nor a character vector")
      else{
        if(any(!is.element(data$region, names(data$population))))
          stop("Some regions are not in the population vector")
        data$population <- 
          c(data$population[c(unique(data$region))],
            data$population[!is.element(names(data$population),
                                        data$region)])

        
        data$distance <- data$distance[names(data$population),
                                       names(data$population)]
        data$region <- as.numeric(factor(data$region, 
                                           levels = unique(data$region)))
      }
    } else{      
      if(length(data$population)<max(data$region))
        stop("The length of the population vector is lower than the maximum value of the region vector")
      if(any(dim(data$distance)<max(data$region)))
        stop("The dimension of the distance matrix is lower than the maximum value of the region vector")
    }
    # can_be_ances_reg is now deprecated
    data$can_be_ances_reg <- matrix(-1, nrow = length(unique(data$region)), 
                                    ncol = length(unique(data$region)))
  }
  if (!is.null(data$s_dens)){
    if(!is.matrix(data$s_dens)) stop("s_dens should be a matrix")
    if(dim(data$s_dens)[1] != dim(data$s_dens)[2]) stop("s_dens should be a square matrix")
    if(dim(data$s_dens)[1] < max(data$region)) stop("Some regions are not in the distance matrix")
    if(any(data$s_dens < 0)) stop("all values in s_dens should be positive")
    if(any(colSums(data$s_dens) != 1)) 
      data$s_dens <- t(t(data$s_dens)/ colSums(data$s_dens))
    if (any(!is.finite(data$s_dens))) {
      stop("non-finite values detected in s_dens")
    }
    
    data$s_dens <- data$s_dens[names(data$population),
                               names(data$population)]
    data$log_s_dens <- (data$s_dens)
  }
  
  ## output is a list of checked data
  return(data)
  
}

