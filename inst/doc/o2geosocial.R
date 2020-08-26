## ----import_library, warning = FALSE, message = FALSE-------------------------
library(o2geosocial)
## We used the data.table syntax throughout this example
library(data.table)
data("toy_outbreak_short")
print(toy_outbreak_short$cases, nrows = 10)
# Extract dataset
dt_cases <- toy_outbreak_short[["cases"]]
dt_cases <- dt_cases[order(Date), ]

## ----plot_ref, warning = FALSE, fig.width = 6.5, fig.height = 4, fig.cap="Figure 2: Cluster size distribution of the simulated dataset"----
# Reference cluster size distribution
hist(table(dt_cases$cluster), breaks = 0:max(table(dt_cases$cluster)), 
     ylab = "Number of clusters", xlab = "Cluster size", main = "", las=1)

## ----distrib------------------------------------------------------------------
# Distribution of the latent period
f_dens <- dgamma(x = 1:100, scale = 0.43, shape = 26.83)
# Distribution of the generation time
w_dens <- dnorm(x = 1:100, mean = 11.7, sd = 2.0)

## ----age, message = FALSE, warning = FALSE------------------------------------
# From the list toy_outbreak_short  
a_dens <- toy_outbreak_short$age_contact
# Alternatively, from POLYMOD:
polymod_matrix <-
  t(socialmixr::contact_matrix(socialmixr::polymod, 
                               countries="United Kingdom",
                               age.limits=seq(0, 70, by=5),
                               missing.contact.age="remove",
                               estimated.contact.age="mean", 
                               missing.participant.age="remove")$matrix)
polymod_matrix <-data.table::as.data.table(polymod_matrix)
# Compute the proportion of connection to each age group
a_dens <- t(t(polymod_matrix)/colSums(polymod_matrix))

## ----distance_pop, warning = FALSE--------------------------------------------
# Extract all regions in the territory
dt_regions <- toy_outbreak_short[["dt_regions"]]
# Extract the population vector
pop_vect <- dt_regions$population
# Create the matrices of coordinates for each region (one "from"; one "to")
mat_dist_from <- matrix(c(rep(dt_regions$long, nrow(dt_regions)),
                          rep(dt_regions$lat, nrow(dt_regions))), ncol = 2)
mat_dist_to <- matrix(c(rep(dt_regions$long, each = nrow(dt_regions)), 
                        rep(dt_regions$lat, each = nrow(dt_regions))),
                      ncol = 2)
# Compute all the distances between the two matrices
all_dist <- geosphere::distGeo(mat_dist_from, mat_dist_to)
# Compile into a distance matrix
dist_mat <- matrix(all_dist/1000, nrow = nrow(dt_regions))
# Rename the matrix columns and rows, and the population vector
names(pop_vect) <- rownames(dist_mat) <- colnames(dist_mat) <- 
  dt_regions$region

## ----first_model, message = FALSE, warning = FALSE----------------------------
# Default moves, likelihoods and priors
moves <- custom_moves()
likelihoods <- custom_likelihoods()
priors <- custom_priors()
# Data and config, model 1
data1 <- outbreaker_data(dates = dt_cases$Date, #Onset dates
                         age_group = dt_cases$age_group, #Age group
                         region = dt_cases$Cens_tract, #Location
                         genotype = dt_cases$Genotype, #Genotype
                         w_dens = w_dens, #Serial interval
                         f_dens = f_dens, #Latent period
                         a_dens = a_dens, #Age stratified contact matrix
                         population = pop_vect, #Population 
                         distance = dist_mat #Distance matrix
)
config1 <- create_config(data = data1, 
                         n_iter = 5000, #Iteration number: main run
                         n_iter_import = 1500, #Iteration number: short run
                         burnin = 500, #burnin period: first run
                         outlier_relative = T, #Absolute / relative threshold 
                         outlier_threshold = 0.9 #Value of the threshold
)
# Run model 1
out1 <- outbreaker(data = data1, config = config1, moves = moves, 
                   priors = priors, likelihoods = likelihoods)
# Set data and config for model 2
data2 <- outbreaker_data(dates = dt_cases$Date, 
                         age_group = dt_cases$age_group,
                         region = dt_cases$Cens_tract,
                         genotype = dt_cases$Genotype, w_dens = w_dens, 
                         f_dens = f_dens, a_dens = a_dens,
                         population = pop_vect, distance = dist_mat,
                         import = dt_cases$import #Import status of the cases
)
config2 <- create_config(data = data2, 
                         find_import = FALSE, # No inference of import status
                         n_iter = 5000, 
                         sample_every = 50, # 1 in 50 iterations is kept
                         burnin = 500)
# Run model 2
out2 <- outbreaker(data = data2, config = config2, moves = moves, 
                   priors = priors, likelihoods = likelihoods)


## ----params, warning = FALSE--------------------------------------------------
burnin <- config1$burnin
# Summary parameters a and b
#Model 1
print(summary(out1, burnin = burnin)$a)
print(summary(out1, burnin = burnin)$b)
# Model 2
print(summary(out2, burnin = burnin)$a)
print(summary(out2, burnin = burnin)$b)

## ----plot_first, warning = FALSE, fig.width = 6.5, fig.height = 4, fig.cap= "Figure 3: Comparison of inferred cluster size distribution with the reference data"----
# We create groups of cluster size: initialise the breaks for each group
group_cluster <- c(1, 2, 3, 5, 10, 15, 40, 100) - 1
# Reference data: h$counts
h <- hist(table(dt_cases$cluster), breaks = group_cluster, plot = FALSE)

# Grouped cluster size distribution in each run
clust_infer1 <- summary(out1, burnin = burnin, 
                        group_cluster = group_cluster)$cluster
clust_infer2 <- summary(out2, group_cluster = group_cluster, 
                        burnin = burnin)$cluster
# Merge inferred and reference cluster size distributions into one matrix
clust_size_matrix <- rbind(clust_infer1["Median",], clust_infer2["Median",],
                           h$counts)
# Plot  
b <- barplot(clust_size_matrix, names.arg = colnames(clust_infer1), las=1,
             ylab = "Number of clusters", xlab = "Cluster size", main = "", 
             beside = T, ylim = c(0, max(c(clust_infer1, clust_infer2))))
# Add the 50% CI
arrows(b[1,], clust_infer1["1st Qu.",], b[1,], clust_infer1["3rd Qu.",], 
       angle = 90, code = 3, length = 0.1)
arrows(b[2,], clust_infer2["1st Qu.",], b[2,], clust_infer2["3rd Qu.",], 
       angle = 90, code = 3, length = 0.1)
# Add legend
legend("topright", fill = grey.colors(3), bty = "n",
       legend = c("Inferred import status", 
                  "Known import status", "Simulated dataset"))

## ----index_infer, warning = FALSE---------------------------------------------
#' Title: Compute the proportion of iterations in the outbreaker() output 
#` where the inferred index matches the actual index in dt_cases
#'
#' @param dt_cases: reference dataset
#' @param out: Matrix output of outbreaker()
#' @param burnin: Numeric, length of the burnin phase
#'
#' @return Numeric vector showing the proportion of iterations pointing to
#' the correct index case
index_infer <- function(dt_cases, out, burnin){
  out_index <- out[out$step > burnin, grep("alpha", colnames(out))]
  ID_index <- matrix(dt_cases[unlist(out_index), ID], ncol = nrow(dt_cases))
  match_infer_data <- t(ID_index) == dt_cases$infector_ID
  match_infer_data[is.na(t(ID_index)) & is.na(dt_cases$infector_ID)] <- TRUE
  prop_correct <- rowSums(match_infer_data, na.rm = T)/ncol(match_infer_data)
  
  return(prop_correct)
}
# Same as index_infer, except it returns the proportion of inferred indexes
# who are in the same reference cluster as the case
index_clust <- function(dt_cases, out, burnin){
  out_index <- out[out$step > burnin, grep("alpha", colnames(out))]
  clust_index <- matrix(dt_cases[unlist(out_index), cluster], 
                        ncol = nrow(dt_cases))
  match_infer_data <- t(clust_index) == dt_cases$cluster
  match_infer_data <- match_infer_data[!is.na(dt_cases$infector_ID),]
  
  
  prop_correct <- rowSums(match_infer_data, na.rm = T)/ncol(match_infer_data)
  
  return(prop_correct)
}
# Run index_infer for each model
index_infer1 <- index_infer(dt_cases = dt_cases, out = out1, burnin = burnin)
index_infer2 <- index_infer(dt_cases = dt_cases, out = out2, burnin = burnin)
# Run index_clust for each model
index_clust1 <- index_clust(dt_cases = dt_cases, out = out1, burnin = burnin)
index_clust2 <- index_clust(dt_cases = dt_cases, out = out2, burnin = burnin)

## ----plot index_infer, warning = FALSE, fig.width = 7, fig.height = 3, fig.cap = "Figure 4: Panel A: Proportion of iterations with the correct index for each case; Panel B: Proportion of iterations where the index is from the correct cluster"----
# Plot the sorted proportion in each model
oldpar <- par(no.readonly =TRUE)
par(bty = "n", mfrow = c(1, 2), mar = c(5,4,2,0), oma = c(0, 0, 0, 0))
# Panel A
plot(sort(index_infer1), type = "l", ylab = "Proportion", xlab = "Case", 
     main =  "A", las=1, col = grey.colors(3)[1], lwd = 3)
lines(sort(index_infer2), col = grey.colors(3)[2], lwd = 3)

# Panel B
plot(sort(index_clust1), type = "l", xlab = "Case", ylab = "", 
     main =  "B", las=1, col = grey.colors(3)[1], lwd = 3)
lines(sort(index_clust2), col = grey.colors(3)[2], lwd = 3)
legend("bottomright", col = grey.colors(3)[1:2], lwd = 3, bty = "n",
       legend = c("Inferred import status", 
                  "Known import status"))

## ----import_sf_file, message = FALSE, warning = FALSE, fig.width = 6.5, fig.height = 2.5----
library(ggplot2)
# Read the shapefile and create one map for each model
map1 <- tigris::tracts(state = "Ohio", class = "sf", progress_bar = FALSE)
map1$INTPTLON <- as.numeric(map1$INTPTLON)
map1$INTPTLAT <- as.numeric(map1$INTPTLAT)
map2 <- map1
map1$model <- "Model 1"
map2$model <- "Model 2"


## ----n_import_reg_per_case, warning = FALSE-----------------------------------
# Add the proportion of iteration in model 1 where each case is an import
dt_cases[, prop_infer1 := summary(out1, burnin = burnin)$tree$import]
# Add the proportion of iteration in model 2 where each case is an import
dt_cases[, prop_infer2 := summary(out2, burnin = burnin)$tree$import]

## ----n_import_reg_per_reg, warning = FALSE------------------------------------
# Number of import per region in model 1
prop_reg1 <- dt_cases[, .(prop_per_reg = sum(prop_infer1)), 
                      by = Cens_tract]$prop_per_reg
# Number of import per region in model 2
prop_reg2 <- dt_cases[, .(prop_per_reg = sum(prop_infer2)), 
                      by = Cens_tract]$prop_per_reg
names(prop_reg1) <- names(prop_reg2) <- unique(dt_cases$Cens_tract)

# Add the number of imports in each region to the maps
map1$prop_reg <- prop_reg1[as.character(map1$GEOID)]
map2$prop_reg <- prop_reg2[as.character(map2$GEOID)]

## ----create_map1, warning = FALSE, fig.width = 6.5, fig.height = 2.5, fig.cap = "Figure 5: Average number of import cases per census tract"----
# Merge maps
maps <- rbind(map1, map2)

limits_long <- c(-84, -82)
limits_lat <- c(40, 41.5)
maps <- maps[maps$INTPTLON > limits_long[1] & maps$INTPTLON < limits_long[2],]
maps <- maps[maps$INTPTLAT > limits_lat[1] & maps$INTPTLAT < limits_lat[2],]

# Plot: number of imports per region, two panels
ggplot(maps) +  geom_sf(aes(fill = prop_reg)) + facet_grid(~model)  +     
  scale_fill_gradient2(na.value = "lightgrey", midpoint = 0.8, 
                       breaks = c(0, 0.5, 1, 1.5), name = "Nb imports",
                       low = "white", mid = "lightblue", high = "darkblue") + 
  coord_sf(xlim = c(-83.8, -82.2), ylim = c(40.2, 41.3)) +
  theme_classic(base_size = 9)

## ----n_sec_region, warning = FALSE--------------------------------------------
#' Title: Compute the number of secondary cases per case in each region
#'
#' @param dt_cases: reference dataset
#' @param out: Matrix output of outbreaker()
#' @param burnin: Numeric, length of the burnin phase
#'
#' @return A numeric matrix: the first column is the census tract ID, the
#' other columns show the number of secondary cases per case. Each row 
#' corresponds to a different iteration.
n_sec_per_reg <- function(dt_cases, out, burnin){
  ## Number of secondary cases per case
  n_sec <- apply(out[out$step > burnin, grep("alpha", colnames(out))], 1, 
                 function(X){
                   X <- factor(X, 1:length(X))
                   return(table(X))})
  ## Aggregate by region
  tot_n_sec_reg <- aggregate(n_sec, list(dt_cases$Cens_tract), sum)
  ## Divide by the number of cases in each region
  tot_n_sec_reg <- cbind(tot_n_sec_reg[, 1], 
                         tot_n_sec_reg[, -1] / table(dt_cases$Cens_tract))
  return(tot_n_sec_reg)
}
  
n_sec_tot1 <- n_sec_per_reg(dt_cases = dt_cases, out = out1, burnin = burnin)
n_sec_tot2 <- n_sec_per_reg(dt_cases = dt_cases, out = out2, burnin = burnin)
n_sec1 <- apply(n_sec_tot1[,-1], 1, median)
n_sec2 <- apply(n_sec_tot2[,-1], 1, median)
names(n_sec1) <- names(n_sec2) <- unique(dt_cases$Cens_tract)

## ----add_variables_sec, warning = FALSE---------------------------------------
map1$n_sec <- as.numeric(n_sec1[as.character(map1$GEOID)])
map2$n_sec <- as.numeric(n_sec2[as.character(map2$GEOID)])

## ----create_maps2, warning = FALSE, fig.width = 6.5, fig.height = 2.5, fig.cap = "Figure 6: Median number of secondary transmission per case in each census tract"----
# Merge maps
maps <- rbind(map1, map2)

limits_long <- c(-84, -82)
limits_lat <- c(40, 41.5)
maps <- maps[maps$INTPTLON > limits_long[1] & maps$INTPTLON < limits_long[2],]
maps <- maps[maps$INTPTLAT > limits_lat[1] & maps$INTPTLAT < limits_lat[2],]

# Plot the geographical distribution of the number of secondary cases
ggplot(maps) +  geom_sf(aes(fill = n_sec)) + facet_grid(~model)  +     
  scale_fill_gradient2(na.value = "lightgrey", mid = "lightblue",
                       low = "white", midpoint = 1, high = "darkblue",
                       breaks = seq(0, 5, 0.5),name = "Sec cases") +
  coord_sf(xlim = c(-83.8, -82.2), ylim = c(40.2, 41.3)) +
  theme_classic(base_size = 9) 

## ----create_stouff_matrix, message = FALSE, warning = FALSE-------------------
# For every column of the distance matrix, use the cumulative sum of the 
# population vector ordered by the distance. Remove the values where 
# the distance between the regions is above gamma
dist_mat_stouffer <- apply(dist_mat, 2, function(X){
  pop_X <- cumsum(pop_vect[order(X)])
  omega_X <- pop_X[names(X)]
  # omega_X is set to -1 if the distance between two regions is above gamma
  omega_X[X > config1$gamma] <- -1
  return(omega_X)
})
# The new value of gamma is equal to the maximum of dist_mat_stouffer + 1
gamma <- max(dist_mat_stouffer) + 1
# The values previously set to -1 are now set to the new value of gamma
dist_mat_stouffer[dist_mat_stouffer == -1] <- max(dist_mat_stouffer) * 2

## ----moves_priors, message = FALSE, warning = FALSE---------------------------
moves3 <- custom_moves(a = cpp_stouffer_move_a)

f_null <- function(param) {
  return(0.0)
}
priors3 <- custom_priors(b = f_null)

## ----out_stouffer, message = FALSE, warning = FALSE---------------------------
# Set data and config lists
data3 <- outbreaker_data(dates = dt_cases$Date, #Onset dates
                         age_group = dt_cases$age_group, #Age group
                         region = dt_cases$Cens_tract, #Location
                         genotype = dt_cases$Genotype, #Genotype
                         w_dens = w_dens, #Serial interval
                         f_dens = f_dens, #Latent period
                         a_dens = a_dens, #Age stratified contact matrix
                         population = pop_vect, #Population 
                         distance = dist_mat_stouffer #Distance matrix
)
config3 <- create_config(data = data3, 
                         gamma = gamma,
                         move_b = FALSE, # b is not estimated
                         init_b = 0, 
                         spatial_method = "power-law",
                         n_iter = 5000, #Iteration number: main run
                         n_iter_import = 1500, #Iteration number: short run
                         burnin = 500, #burnin period: first run
                         outlier_relative = T, #Absolute / relative threshold
                         outlier_threshold = 0.9 #Value of the threshold
)
# Run the model using the Stouffer's rank method
out_stouffer <- outbreaker(data = data3, config = config3, moves = moves3, 
                           priors = priors3, likelihoods = likelihoods)

## ----stouffer_cluster, message = FALSE, warning = FALSE, fig.width = 6.5, fig.height = 4, fig.cap = "Figure 7: Comparison of inferred cluster size distribution with the reference data"----
# Grouped cluster size distribution in the Stouffer's rank model
clust_infer_stouf <- summary(out_stouffer, burnin = burnin, 
                             group_cluster = group_cluster)$cluster
# Merge inferred and reference cluster size distributions
clust_size_matrix <- rbind(clust_infer_stouf["Median",], h$counts) 
# Plot the two distributions
b <- barplot(clust_size_matrix, names.arg = colnames(clust_infer_stouf), 
             beside = T)
# Add CIs
arrows(b[1,], clust_infer_stouf["1st Qu.",], b[1,], 
       clust_infer_stouf["3rd Qu.",], angle = 90, code = 3, length = 0.1)
legend("topright", fill = grey.colors(2), bty = "n",
       legend = c("Inferred import status, Stouffer's rank method", 
                  "Simulated dataset"))

## ----match_stouffer, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 3, fig.cap = "Figure 8: Panel A: Proportion of iterations with the correct index for each case; Panel B: Proportion of iterations where the index is from the correct cluster"----
index_infer_stouf <- index_infer(dt_cases = dt_cases, out = out_stouffer, 
                                 burnin = config1$burnin)
index_clust_stouf <- index_clust(dt_cases = dt_cases, out = out_stouffer, 
                                 burnin = config1$burnin)
# Plot the sorted proportion in each model
par(bty = "n", mfrow = c(1, 2), mar = c(5,4,2,0), oma = c(0, 0, 0, 0))
# Panel A
plot(sort(index_infer_stouf), type = "l", xlab = "Case", ylab = "", 
     main =  "A", las=1, col = grey.colors(2)[1], lwd = 3)
# Panel B
plot(sort(index_clust_stouf), type = "l", ylab = "Proportion", xlab = "Case", 
     main =  "B", las=1, col = grey.colors(2)[1], lwd = 3)
par(oldpar)

