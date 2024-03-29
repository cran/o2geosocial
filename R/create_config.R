#' Set and check parameter settings 
#'
#' This function defines settings.  It takes a list of named
#' items as input, performs various checks, set defaults where arguments are
#' missing, and return a correct list of settings. If no input is given, it
#' returns the default settings.
#'
#' Acceptable arguments for ... are:
#'
#' \describe{
#'
#' \item{init_tree}{the tree used to initialize the MCMC. It can be a vector of
#' integers corresponding to the tree itself, where the i-th value corresponds
#' to the index of case \code{i}. Otherwise, it should be defined as the character 
#' string "star" and the function \code{create_config()} will generate the initial 
#' tree.}
#'
#' \item{spatial_method}{a character string indicating the method used to
#' evaluate the spatial likelihood. Can be either "exponential" or "power-law".}
#'
#' \item{gamma}{a double indicating the spatial threshold for pre clustering; 
#' defaults to NULL.}
#' 
#' \item{delta}{a double indicating the temporal threshold for pre clustering;  
#' defaults to NULL.}
#' 
#' \item{init_kappa}{a vector of integers indicating the initial values of kappa; 
#' defaults to 1.}
#' 
#' \item{init_a}{initial value of the first spatial parameter (population).}
#' 
#' \item{init_b}{initial value of the second spatial parameter (distance).}
#' 
#' \item{init_alpha}{a vector of integers indicating the initial values of
#' alpha, where the i-th value indicates the ancestor of case 'i'; defaults to
#' \code{NULL}, in which ancestries are defined from \code{init_tree}.}
#'
#' \item{init_t_inf}{a vector of integers indicating the initial values of
#' \code{t_inf}, i.e. dates of infection; defaults to \code{NULL}, in which case
#' the most likely \code{t_inf} will be determined from the delay to
#' reporting/symptoms distribution, and the dates of reporting/symptoms,
#' provided in \code{data}.}
#'
#' \item{init_pi}{initial value for the reporting probability.}
#'
#' \item{n_iter}{an integer indicating the number of iterations in the MCMC,
#' including the burnin period.}
#'
#' \item{move_alpha}{a vector of logicals indicating, for each case, if the
#' ancestry should be estimated ('moved' in the MCMC), or not, defaulting to
#' TRUE; the vector is recycled if needed.}
#'
#' \item{move_t_inf}{a vector of logicals indicating, for each case, if the
#' dates of infection should be estimated ('moved' in the MCMC), or not,
#' defaulting to TRUE; the vector is recycled if needed.}
#'
#' \item{move_pi}{a logical indicating whether the reporting probability
#' should be estimated ('moved' in the MCMC), or not, all defaulting to TRUE.}
#'
#' \item{move_kappa}{a logical indicating whether the number of generations
#' between two successive cases should be estimated ('moved' in the MCMC), or
#' not, all defaulting to TRUE.}

#' \item{move_a}{a logical indicating whether the first spatial parameter
#' should be estimated ('moved' in the MCMC), or not, all defaulting to TRUE.}

#' \item{move_b}{a logical indicating whether the second spatial parameter
#' should be estimated ('moved' in the MCMC), or not, all defaulting to TRUE.}
#' 
#' \item{move_swap_cases}{a logical indicating whether the movement to swap cases
#' should be used, or not, all defaulting to TRUE.}
#'
#' \item{sample_every}{the frequency at which MCMC samples are retained for the
#' output.}
#'
#' \item{sd_pi}{the standard deviation for the Normal proposal for the reporting
#' probability.}
#' 
#' \item{sd_a}{the standard deviation for the Normal proposal for the first spatial 
#' parameter.}
#' 
#' \item{sd_b}{the standard deviation for the Normal proposal for the second spatial 
#' parameter.}
#'
#' \item{find_import}{a logical indicating whether the import status of cases should
#' be estimated.}
#' 
#' \item{outlier_threshold}{a numeric value indicating the probability that should be
#' used to compute the threshold when estimating the import status.}
#'
#' \item{outlier_relative}{a logical indicating whether the threshold is an absolute 
#' or relative value, default to FALSE (absolute value).}
#' 
#' \item{outlier_plot}{a logical indicating whether to plot the comparison between the 
#' likelihoods of connection in the short run and the threshold.}
#' 
#' \item{n_iter_import}{Number of iterations of the first short run.}
#' 
#' \item{sample_every_import}{the frequency at which MCMC samples are retained for the
#' output during the first run.}
#' 
#' \item{burnin}{The number of iterations that should be removed when estimating import.}
#' 
#' \item{max_kappa}{an integer indicating the largest number of generations
#' between any two linked cases; defaults to 5.}
#'
#' \item{prior_pi}{a numeric vector of length 2 indicating the first and second
#' parameter of the beta prior for the reporting probability 'pi'.}
#' 
#' \item{prior_a}{a numeric vector of length 2 indicating the first and second
#' parameter of the uniform prior for the first spatial parameter 'a'.}
#' 
#' \item{prior_b}{a numeric vector of length 2 indicating the first and second
#' parameter of the uniform prior for the second spatial parameter 'b'.}

#' \item{verbatim}{Logical, should the number of iteration be printed.}
#'
#' }
#'
#' @param data an optional list of data items as returned by
#'     \code{outbreaker_data}; if provided, this allows for further checks of
#'     the outbreaker settings.
#'     
#' @param ... a list of config items to be processed (see description)
#' 
#' @seealso \code{\link{outbreaker_data}} to check and process data for outbreaker
#'
#' @author Initial version by Thibaut Jombart, rewritten by Alexis Robert (\email{alexis.robert@lshtm.ac.uk})
#'
#' @return 
#' 
#' A named list containing the value of each elements listed in the 'Details' section. This list describes the settings of the \code{outbreaker()} function. The class of this list is set to \code{outbreaker_config}.
#' @export
#'
#' @examples
#' ## see default settings
#' create_config()
#'
#' ## change defaults
#' create_config(move_alpha = FALSE, n_iter = 2e5, sample_every = 1000)
#'
#'
#'
create_config <- function (..., data = NULL) 
{
  config <- list(...)
  if (length(config) == 1L && is.list(config[[1]])) {
    config <- config[[1]]
  }
  defaults <- list(init_tree = c("star"), spatial_method = "exponential",
                   gamma = NULL,delta = NULL,
                   init_alpha = NULL, init_kappa = 1, init_t_inf = NULL, 
                   init_pi = 0.9, init_a = 1, init_b = 0.1,
                   
                   move_alpha = TRUE, move_swap_cases = TRUE, 
                   move_t_inf = TRUE, move_kappa = TRUE, move_pi = TRUE, 
                   move_a = TRUE, move_b = TRUE,
                   
                   prior_pi = c(10, 1), 
                   prior_a = c(0.2, 1.5), prior_b = c(0.01, 2), 
                   
                   sd_pi = 0.1, sd_a = 0.1, sd_b = 0.1, 
                   
                   n_iter = 10000, sample_every = 50, 
                   
                   max_kappa = 2, find_import = TRUE, outlier_threshold = 0.05, 
                   outlier_relative = FALSE, outlier_plot = FALSE,
                   
                   n_iter_import = 5000, sample_every_import = 50, 
                   burnin = 1000, verbatim = FALSE)
  config <- modify_defaults(defaults, config)
  if (is.character(config$init_tree)) {
    config$init_tree <- match.arg(config$init_tree, c("star"))
  }
  if (config$spatial_method == "exponential" ||
      config$spatial_method == "power-law") {
    config$spatial_method <- match.arg(config$spatial_method, c("exponential", 
                                                                "power-law"))
  } else if(!is.null(config$spatial_method) || (is.null(config$spatial_method) &&
                                                (config$move_a == TRUE || 
                                                 config$move_b == TRUE))){
    stop("invalid value for spatial_method, spatial_method is either exponential, or power-law.")
  } 
  
  if(!is.null(config$gamma)){
    if(is.na(config$gamma))
      stop("gamma is NA")
    if(config$gamma<= 0)
      stop("gamma is below 0")
    if(!is.numeric(config$gamma))
      stop("gamma is not numeric")
  }
  if(config$gamma == 0 || is.null(config$gamma)){
    config$gamma <- Inf
  }
  if(!is.null(config$delta)){
    if(is.na(config$delta))
      stop("delta is NA")
    if(config$delta <= 0)
      stop("delta is below 0")
    if(!is.numeric(config$delta))
      stop("delta is not numeric")
  }
  if(config$delta == 0 || is.null(config$delta)){
    config$delta <- Inf
  }
  if (is.numeric(config$init_tree)) {
    tree_ances <- config$init_tree
    if(!any(is.na(tree_ances))) 
      stop("There should be an ancestor in the initial tree")
    i <- 0
    while(any(!is.na(tree_ances[tree_ances]))){ 
      tree_ances[!is.na(tree_ances[tree_ances])] <- 
        tree_ances[tree_ances][!is.na(tree_ances[tree_ances])]
      i <- i + 1
      if(i > 20) stop("invalid initial tree")
    }
    tree_ances[is.na(tree_ances)] <- which(is.na(tree_ances))
    genotype_tree <- numeric(length(unique(tree_ances)))
    nb_gen_rep_per_tree <- vapply(unique(tree_ances), function(X) {
      gens <- unique(data$genotype[which(tree_ances == X)])
      return(length(gens[gens != "Not attributed"]))
    }, 1)
    if(any(nb_gen_rep_per_tree >= 2))
      stop("There should be one reported genotype per tree at most.")
    config$init_alpha <- as.integer(config$init_tree)
  }
  if (!is.null(config$init_t_inf)) {
    if (inherits(config$init_t_inf, "Date")) {
      config$init_t_inf <- config$init_t_inf - min(config$init_t_inf)
    }
    if (inherits(config$init_t_inf, "POSIXct")) {
      config$init_t_inf <- difftime(config$init_t_inf, 
                                    min(config$init_t_inf), units = "days")
    }
    config$init_t_inf <- as.integer(round(config$init_t_inf))
  }
  if (!is.null(config$init_alpha)) {
    are_not_imports <- !is.na(config$init_alpha)
  }
  else {
    are_not_imports <- TRUE
  }
  if (!is.numeric(config$init_kappa)) {
    stop("init_kappa is not a numeric value")
  }
  config$init_kappa <- as.integer(round(config$init_kappa))
  if (any(config$init_kappa < 1, na.rm = TRUE)) {
    stop("init_kappa has values smaller than 1")
  }
  if (any(config$init_kappa > config$max_kappa, na.rm = TRUE)) {
    config$init_kappa[config$init_kappa > config$max_kappa] <- config$max_kappa
    warning("values of init_kappa greater than max_kappa have been set to max_kappa")
  }
  if (!is.numeric(config$init_pi)) {
    stop("init_pi is not a numeric value")
  }
  if (!is.finite(config$init_pi)) {
    stop("init_pi is infinite or NA")
  }
  if (config$init_pi < 0) {
    stop("init_pi is negative")
  }
  if (config$init_pi > 1) {
    stop("init_pi is greater than 1")
  }
  if(!is.numeric(config$init_a)){
    stop("init_a is not a numeric value")
  }
  if (config$init_a < 0) {
    stop("init_a is negative")
  }
  if (!is.finite(config$init_a)) {
    stop("init_a is infinite or NA")
  }
  if(!is.numeric(config$init_b)){
    stop("init_b is not a numeric value")
  }
  if (config$init_b < 0) {
    stop("init_b is negative")
  }
  if (!is.finite(config$init_b)) {
    stop("init_b is infinite or NA")
  }
  if (!all(is.logical(config$move_alpha))) {
    stop("move_alpha is not a logical")
  }
  if (any(is.na(config$move_alpha))) {
    stop("move_alpha is NA")
  }
  if (!is.logical(config$move_swap_cases)) {
    stop("move_swap_cases is not a logical")
  }
  if (is.na(config$move_swap_cases)) {
    stop("move_swap_cases is NA")
  }
  if (!is.logical(config$move_t_inf)) {
    stop("move_t_inf is not a logical")
  }
  if (any(is.na(config$move_t_inf))) {
    stop("move_t_inf has NAs")
  }
  if (!is.logical(config$move_kappa)) {
    stop("move_kappa is not a logical")
  }
  if (any(is.na(config$move_kappa))) {
    stop("move_kappa has NA")
  }
  if (!is.logical(config$move_pi)) {
    stop("move_pi is not a logical")
  }
  if (is.na(config$move_pi)) {
    stop("move_pi is NA")
  }
  if (!is.logical(config$move_a)) {
    stop("move_a is not a logical")
  }
  if (is.na(config$move_a)) {
    stop("move_a is NA")
  }
  if (!is.logical(config$move_b)) {
    stop("move_b is not a logical")
  }
  if (is.na(config$move_b)) {
    stop("move_b is NA")
  }
  if (!is.numeric(config$n_iter)) {
    stop("n_iter is not a numeric value")
  }
  if (!is.numeric(config$burnin)) {
    stop("burnin is not a numeric value")
  }
  if (config$n_iter < 2) {
    stop("n_iter is smaller than 2")
  }
  if (!is.finite(config$n_iter)) {
    stop("n_iter is infinite or NA")
  }
  if (config$burnin < 0) {
    stop("n_iter is negative")
  }
  if (!is.finite(config$burnin)) {
    stop("burnin is infinite or NA")
  }
  if (!is.numeric(config$sample_every)) {
    stop("sample_every is not a numeric value")
  }
  if (config$sample_every < 1) {
    stop("sample_every is smaller than 1")
  }
  if (!is.finite(config$sample_every)) {
    stop("sample_every is infinite or NA")
  }
  config$sample_every <- as.integer(config$sample_every)
  if (!is.numeric(config$sd_pi)) {
    stop("sd_pi is not a numeric value")
  }
  if (config$sd_pi < 1e-10) {
    stop("sd_pi is close to zero or negative")
  }
  if (!is.finite(config$sd_pi)) {
    stop("sd_pi is infinite or NA")
  }
  if (!is.numeric(config$sd_a)) {
    stop("sd_a is not a numeric value")
  }
  if (config$sd_a < 1e-10) {
    stop("sd_a is close to zero or negative")
  }
  if (!is.finite(config$sd_a)) {
    stop("sd_a is infinite or NA")
  }
  if (!is.numeric(config$sd_b)) {
    stop("sd_b is not a numeric value")
  }
  if (config$sd_b < 1e-10) {
    stop("sd_b is close to zero or negative")
  }
  if (!is.finite(config$sd_b)) {
    stop("sd_b is infinite or NA")
  }
  if (!is.logical(config$find_import)) {
    stop("find_import is not logical")
  }
  if (length(config$find_import) != 1L) {
    stop("find_import should be a single logical value")
  }
  if (is.na(config$find_import)) {
    stop("find_import is NA")
  }
  if (!is.numeric(config$outlier_threshold)) {
    stop("outlier_threshold is not a numeric value")
  }
  if (!is.logical(config$outlier_relative)) {
    stop("outlier_relative is not a logical value")
  }
  if (!is.finite(config$outlier_threshold)) {
    stop("outlier_threshold is infinite or NA")
  }
  if (!is.logical(config$outlier_plot)) {
    stop("outlier_plot is not a logical")
  }
  if (is.na(config$outlier_plot)) {
    stop("outlier_plot is NA")
  }
  if (!is.numeric(config$n_iter_import)) {
    stop("n_iter_import is not a numeric value")
  }
  if (!is.finite(config$n_iter_import)) {
    stop("n_iter_import is infinite or NA")
  }
  config$n_iter_import <- as.integer(config$n_iter_import)
  if (!is.numeric(config$sample_every_import)) {
    stop("sample_every_import is not a numeric value")
  }
  if (config$sample_every_import < 1) {
    stop("sample_every_import is smaller than 1")
  }
  if (!is.finite(config$sample_every_import)) {
    stop("sample_every_import is infinite or NA")
  }
  config$sample_every_import <- as.integer(config$sample_every_import)
  if (!all(is.numeric(config$prior_pi))) {
    stop("prior_pi has non-numeric values")
  }
  if (any(config$prior_pi < 0)) {
    stop("prior_pi has negative values")
  }
  if (length(config$prior_pi) != 2L) {
    stop("prior_pi should be a vector of length 2")
  }
  if (!all(is.finite(config$prior_pi))) {
    stop("prior_pi is has values which are infinite or NA")
  }
  if (!all(is.numeric(config$prior_a))) {
    stop("prior_a has non-numeric values")
  }
  if (any(config$prior_a < 0)) {
    stop("prior_a has negative values")
  }
  if (length(config$prior_a) != 2L) {
    stop("prior_a should be a vector of length 2")
  }
  if (!all(is.finite(config$prior_a))) {
    stop("prior_a is has values which are infinite or NA")
  }
  if (!all(is.numeric(config$prior_b))) {
    stop("prior_b has non-numeric values")
  }
  if (any(config$prior_b < 0)) {
    stop("prior_b has negative values")
  }
  if (length(config$prior_b) != 2L) {
    stop("prior_b should be a vector of length 2")
  }
  if (!all(is.finite(config$prior_b))) {
    stop("prior_b is has values which are infinite or NA")
  }
  if(!is.null(config$verbatim)){
    if(!is.logical(config$verbatim)) stop("verbatim should be a logical value")
  }
  
  if((config$move_a == TRUE || config$move_b == TRUE) && config$max_kappa > 2){
    warning("If spatial kernel parameters are estimated, max_kappa is set to 2")
    config$max_kappa <- 2
  }
  if (!is.null(data)) {
    if (is.character(config$init_tree)) {
      ## Config$init_tree is not a vector of numeric values (showing the initial trees),
      ## the function has to generate the initial trees with max 1 genotype reported per tree
      config$init_alpha <- data$import*0
      config$init_alpha[data$import == TRUE] <- NA
      config$init_alpha[!duplicated(data$is_cluster)] <- NA
      for(X in unique(data$is_cluster)){
        # Cases in the cluster
        cases_clust <- which(data$is_cluster == X)
        # All genotypes in the cluster
        nb_gen <- unique(data$genotype[data$is_cluster == X & 
                                         data$genotype != "Not attributed"])
        # Genotypes of cases that are imports
        gen_clust <- data$genotype[cases_clust][is.na(config$init_alpha[cases_clust])]
        # Unique genotypes + all NAs
        unique_gen_clust <- gen_clust[gen_clust == "Not attributed" | 
                                        !duplicated(gen_clust)]
        # If more genotypes than potential imports, add new import
        while(length(nb_gen) > length(unique_gen_clust)){
          config$init_alpha[data$is_cluster == X & 
                              data$genotype != "Not attributed" &
                              !is.element(data$genotype, unique_gen_clust)][1] <- NA
          gen_clust <- data$genotype[cases_clust][is.na(config$init_alpha[cases_clust])]
          unique_gen_clust <- gen_clust[gen_clust == "Not attributed" | 
                                          !duplicated(gen_clust)]
          
        }
        # All genotypes in non imports
        nb_gen_sec <- unique(data$genotype[data$is_cluster == X & 
                                             !is.na(config$init_alpha) &
                                             data$genotype != "Not attributed"])
        pot_ances <- which(is.na(config$init_alpha) & data$is_cluster == X)
        count <- 1
        for(j in nb_gen_sec){
          count_gen <- count
          notances_X_j <- which(data$is_cluster == X & 
                                  !is.na(config$init_alpha) &
                                  data$genotype == j)
          
          # If there is an ancestor from X with the same genotype before the first case
          # then increment count_gen
          if(any(is.na(config$init_alpha) & data$is_cluster == X & 
                 data$genotype == j & 
                 data$dates <= min(data$dates[notances_X_j]))){
            count_gen <- which(gen_clust == j)[1]
          } else {
            while(count_gen < length(gen_clust) && !(gen_clust[count_gen] == j)){
              count_gen <- count_gen + 1
            }
            if(gen_clust[count_gen] != j){
              count_gen <- count
              while(count_gen < length(gen_clust) &&
                    !(gen_clust[count_gen] == "Not attributed")){
                count_gen <- count_gen + 1
              }
            }
            if(gen_clust[count_gen] == "Not attributed") count <- count_gen + 1
          }
          # If the potential import was reported before all the cases, then the first
          # case is added as an import
          if(data$dates[pot_ances[count_gen]] > min(data$dates[notances_X_j])){
            new_ances <- notances_X_j[1]
            config$init_alpha[new_ances] <- NA
            config$init_alpha[data$is_cluster == X & !is.na(config$init_alpha) &
                                data$genotype == j] <- new_ances
          } else {
            # Otherwise, all cases are linked to the potential import
            config$init_alpha[notances_X_j] <- pot_ances[count_gen]
          }
        }
        for (k in which(data$is_cluster == X &
                        !is.na(config$init_alpha) &
                        data$genotype == "Not attributed")){
          config$init_alpha[k] <- max(which(data$is_cluster[1:k] == X &
                                              is.na(config$init_alpha[1:k])))
        }
      }
      
    } else {
      if (length(config$init_alpha) != data$N) {
        stop("length of init_alpha or init_tree incorrect")
      }
      unknownAnces <- config$init_alpha < 1 | config$init_alpha > 
        data$N
      if (any(stats::na.omit(unknownAnces))) {
        warning("some initial ancestries refer to unknown cases (idx<1 or >N)")
        config$init_alpha[unknownAnces] <- NA
      }
    }
    if (!is.null(config$init_t_inf) && !is.null(data$dates)) {
      if (any(config$init_t_inf >= data$dates, na.rm = TRUE)) {
        msg <- paste0("Initial dates of infection come after ", 
                      "sampling dates / dates of onset.")
        stop(msg)
      }
    } 
    else if(!is.null(data$dates)){
      max_like_delay <- which.max(data$f_dens)
      if (!is.finite(max_like_delay) || length(max_like_delay) == 0) {
        max_like_delay <- 1L
      }
      config$init_t_inf <- as.integer(data$dates - max_like_delay)
    }
    else
      config$init_t_inf <- NULL
    config$move_alpha <- rep(config$move_alpha, length.out = data$N)
    config$move_t_inf <- rep(config$move_t_inf, length.out = data$N)
    config$move_kappa <- rep(config$move_kappa, length.out = data$N)
    config$init_kappa <- rep(config$init_kappa, length.out = data$N)
    config$init_kappa[is.na(config$init_alpha)] <- NA
    
    config$move_alpha[data$import == TRUE] <- FALSE
    config$init_kappa[data$import == TRUE] <- NA
    config$init_t_inf[is.na(config$init_alpha)] <-
      config$init_t_inf[is.na(config$init_alpha)]-1
    
    config$init_alpha <- as.integer(config$init_alpha)
    config$init_t_inf <- as.integer(config$init_t_inf)
    if(is.null(data$population) || is.null(data$region) ||
       is.null(data$distance)){
      warning("move_a and move_b are set to FALSE since data$population, data$region or data$distance is null")
      config$move_a <- FALSE
      config$move_b <- FALSE
      config$gamma <- Inf
    }
    if(all(is.na(config$init_alpha)))
      warning("All cases are currently set as imports")
  }
  class(config) <- c("outbreaker_config", "list")
  return(config)
}
