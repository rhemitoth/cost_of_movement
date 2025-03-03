## ---------------------------
##
## Script name: COM test
##
## Author: Rhemi Toth
##
## Date Created: 2025-02-26
##
## Email: rhemitoth@g.harvard.edu
##
## ---------------------------
##
## Notes: testing out using TGP + ARC model to calculate movement costs
## for a roe deer in Cembra
##   
##
## ---------------------------

rm(list = ls())


# Packages ----------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(tgp)
library(coda)
library(terra)
library(sf)

# Load Data ---------------------------------------------------------------

# GPS data
gps_full <- read_csv("/Users/rhemitoth/Library/CloudStorage/GoogleDrive-rhemitoth@g.harvard.edu/My Drive/Cembra/Movement_Data/GPS/gps_data.csv")
dem <- rast("/Users/rhemitoth/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/Rhemi/Cembra/GIS/Images/TINITALY/w51065_s10.tif")

# Prepare the GPS data ----------------------------------------------------

# Select only relevant columns
# Round "hourly" timetamps to the nearest hour to get rid of seconds info

gps <- gps_full %>%
  filter(animals_id == 3540,
         utm_x != "NULL")%>%
  rename(X = utm_x,
         Y = utm_y)%>%
  mutate(timestamp = round_date(acquisition_time, unit = "hour"),
         X = as.numeric(X),
         Y = as.numeric(Y))%>%
  select(timestamp,
         X,
         Y)%>%
  arrange(timestamp)%>%
  filter(date(timestamp) == "2016-12-04")

# Create a t.obs column (hours since first recorded GPS point)
gps <- gps %>%
  mutate(t.obs = as.numeric(timestamp - gps$timestamp[1])/3600) 

# Fit the TGP model -------------------------------------------------------

model_x <- btgpllm(X=gps$t.obs, #fit time as explanatory variable
                   Z= gps$X, #fit easting location as dependent variable
                   bprior="b0", #set uninformative hierarchical normal prior
                   verb=1, #don't print progress meter. Set =1 to print progress meter.
                   pred.n = F) #prevents prediction at t.obs points for faster model

# before updating BLAS: 00:01:40

#fitting (all prediction will be made in later step)

#same as above settings, but fit to northing location:
model_y <- btgpllm(X=gps$t.obs,
                   Z= gps$Y,
                   bprior="b0",
                   verb=1 ,
                   pred.n = F)

# before updating BLAS: 00:04:42


# Predict Paths -----------------------------------------------------------

# Set MCMC hyperparameters
bte <- c(400,500,2) #EDIT this line to prototype
MCMC_sample_size <- (bte[2] - bte[1])/bte[3] #do not edit this line

# Set delta t
delta.t <- 1/12 # 5 minutes

# Use fitted model to predict
season.end <- max(gps$t.obs)
XX <- seq(0, season.end, by = delta.t) #points we want prediction at
#### build storage table
all_paths <- data.frame(matrix(ncol = 1+ 2*MCMC_sample_size, nrow = length(XX)) )
names(all_paths)[1] <- "t"
all_paths$t <- XX
seq.x <- seq(2,ncol(all_paths), by = 2) #sequence of even numbers
names(all_paths)[seq.x] <- paste0('x', 1:MCMC_sample_size)
seq.y <- seq(3,ncol(all_paths), by = 2) #sequence of odd numbers not 1
names(all_paths)[seq.y] <- paste0('y', 1:MCMC_sample_size)
#### tgp model fitting transformation
# The tgp model fitting utilizes a specific transformation of the data used to fit
#the model. This ensures that prespecified priors, the scale hyperparameter, etc.
#set within the package work for all types of data. The specific transformation
#used by the package is to transform data that have a mean of 0 and range of 1.
## make your transformation constants from the data once:
z.x <-gps$X
z.x.ranged <- -0.5 +(z.x - min(z.x))/(max(z.x) - min(z.x))
z.y <-gps$Y
z.y.ranged <- -0.5 +(z.y - min(z.y))/(max(z.y) - min(z.y))
#### compute predictions
predicted_x <- predict(model_x, XX, pred.n=FALSE, BTE=bte, MAP=FALSE, trace=T)

#Before BLAS: 00:23:000

posterior_predictives_x_ <- predicted_x$trace$preds$ZZ
predicted_y <- predict(model_y, XX, pred.n=FALSE, BTE=bte, MAP=FALSE, trace=T)
posterior_predictives_y <- predicted_y$trace$preds$ZZ
#Setting MAP = FALSE when making prediction means that the MCMC sampling will start at
#the MAP (maximum a posteriori probability) tree, but then continue in Bayesian form
#and sample across the full posterior distribution of trees. This means different MCMC
#samples will often be from different possible trees (different partitioning solutions).
# Back transform to undue model fitting scaling:
posterior_predictives_x <- t(posterior_predictives_x_ + mean(z.x.ranged) + 0.5)*
  (max(z.x) - min(z.x)) + min(z.x);
posterior_predictives_y <- t(posterior_predictives_y + mean(z.y.ranged) + 0.5)*
  (max(z.y) - min(z.y)) + min(z.y);
# take transpose of matrix so it's set up with times are rows and samples are columns
#instead of vice versa
### seperate these into the table
all_paths[,seq.x] <- posterior_predictives_x ; #fill in all the xs
all_paths[,seq.y] <- posterior_predictives_y #fill in all the ys

# check the effective sample size
# A large effective sample size (close to the MCMC sample size) implies
# good mixing of the sample and low autocorrelation. Low effective sample size suggests autocorrelation that
# may require increased thinning (increased E in BTE settings).
sample <- as.numeric(posterior_predictives_x_[,1])
effectiveSize(sample)

# trace plot of MCMC samples
# The chain should stabilize around a central value, meaning it has reached the posterior distribution.
# If there is a visible trend (e.g., increasing or decreasing over time), the chain has not converged.
# Good mixing means the chain moves up and down frequently rather than being "sticky" in certain regions.
# Poor mixing suggests high autocorrelation, meaning successive samples are highly dependent on each other.
# If samples are highly correlated, the chain explores the posterior slowly.
# Ideally, the samples should resemble white noise.
# A good trace plot has a fuzzy caterpillar shape
# Bad trace plot: A clear upward or downward trend, strong periodicity, or long stretches where the chain stays in one region.
plot(1:MCMC_sample_size, sample,type="l")


# Compute Mean Location and 95% CI ----------------------------------------

# Convert to long format
all_paths_long <- all_paths %>%
  pivot_longer(-t, names_to = c(".value", "sample"), names_pattern = "([xy])(\\d+)")

# Get mean location and CI
# Compute summary statistics
trajectory <- all_paths_long %>%
  group_by(t) %>%
  summarise(
    x_mean = mean(x, na.rm = TRUE),
    x_95_CI_low = quantile(x, 0.025, na.rm = TRUE),
    x_95_CI_high = quantile(x, 0.975, na.rm = TRUE),
    y_mean = mean(y, na.rm = TRUE),
    y_95_CI_low = quantile(y, 0.025, na.rm = TRUE),
    y_95_CI_high = quantile(y, 0.975, na.rm = TRUE)
  )%>%
  mutate(x_CI = x_95_CI_high - x_mean,
         y_CI = y_95_CI_high - y_mean)


# Plot the x trajectory ---------------------------------------------------

x_traj <- ggplot()+
  geom_ribbon(data = trajectory,
              aes(x = t,
                  ymin =  x_95_CI_low,
                  ymax = x_95_CI_high),
              fill = "blue",
              alpha = 0.4)+
  geom_point(data = trajectory,
             aes(x = t, y = x_mean),
             color = "blue",
             size = 0.5)+
  theme_bw()

x_traj

# Plot the y trajectory ---------------------------------------------------

y_traj <- ggplot()+
  geom_ribbon(data = trajectory,
              aes(x = t,
                  ymin =  y_95_CI_low,
                  ymax = y_95_CI_high),
              fill = "blue",
              alpha = 0.4)+
  geom_point(data = trajectory,
             aes(x = t, y = y_mean),
             color = "blue",
             size = 0.5)+
  theme_bw()

y_traj



# Plot a random sample of trajectories ------------------------------------

# Randomly sample 20 trajectories

# Get trajectory column names (excluding "t")
trajectory_cols <- colnames(all_paths)[-1]

# Extract x-y pairs
trajectory_pairs <- matrix(trajectory_cols, ncol = 2, byrow = TRUE)

# Randomly select 20 trajectories
set.seed(123)  # For reproducibility
selected_indices <- sample(nrow(trajectory_pairs), 20, replace = FALSE)
selected_trajectories <- trajectory_pairs[selected_indices, ]

# Convert selected matrix rows to a vector
selected_trajectories <- as.vector(trajectory_pairs[selected_indices, ])

# Convert data to long format for ggplot
sampled_trajectories <- all_paths %>%
  select(t, all_of(selected_trajectories)) %>%
  pivot_longer(cols = -t,
               names_to = c(".value", "trajectory"),
               names_pattern = "([xy])([0-9]+)")

# Plot

ggplot() +
  geom_path(data=sampled_trajectories,
            aes(x = x, y = y, group = trajectory),
            color = "grey",
            size = 0.2,
            alpha = 0.5) +
  geom_path(data = trajectory,
            aes(x = x_mean, y = y_mean),
            color = "blue",
            size = 0.5,)+
  geom_point(data = gps,
             aes(x = X,
                 y = Y),
             shape = 21,
             fill = "white")+
  labs(title = "Modeled vs Observed Trajectory", x = "X", y = "Y") +
  theme_minimal()


# Plot error elipses around trajectory points -----------------------------

make_ellipse <- function(x, y, x_err, y_err, n = 100) {
  t <- seq(0, 2 * pi, length.out = n)
  data.frame(
    x = x + x_err * cos(t),
    y = y + y_err * sin(t)
  )
}

# Generate ellipses for each point
ellipses <- bind_rows(lapply(1:nrow(trajectory), function(i) {
  make_ellipse(trajectory$x_mean[i], trajectory$y_mean[i], trajectory$x_CI[i], trajectory$y_CI[i]) %>%
    mutate(id = i)
}))

# Plot the trajectory with error ellipses

ggplot() +
  geom_polygon(data = ellipses, aes(x, y, group = id), fill = "blue", alpha = 0.2) +  # Error ellipses
  geom_point(data = trajectory, aes(x_mean, y_mean), fill = "white", size = 2) +  # GPS points
  geom_point(data = gps,
             aes(x = X,
                 y = Y),
             shape = 21,
             fill = "white")+
  theme_minimal() +
  labs(title = "Animal Trajectory with 95% Confidence Ellipses",
       x = "Easting", y = "Northing")



# Compute ECOT using ARC Model --------------------------------------------

# Extract elevation at each point
all_paths_long$z <- extract(dem, all_paths_long[, c("x", "y")])[, 2]

# Convert all_paths_long to a matrix
all_paths_wide <- all_paths_long %>%
  pivot_wider(names_from = sample, values_from = c(x, y, z), names_prefix = "") %>%
  select(t, 
         all_of(
           colnames(.)[-1][order(
             as.numeric(str_extract(colnames(.)[-1], "\\d+")),  # Extract numbers and sort numerically
             match(str_extract(colnames(.)[-1], "[xyz]"), c("x", "y", "z")) # Ensure x, y, z ordering
           )]
         )
  )



compute_distances <- function(all_paths_wide, MCMC_sample_size) {
  n_time <- nrow(all_paths_wide) - 1  # Number of time steps
  
  # Create a correctly indexed data frame
  all_distances <- data.frame(matrix(NA, nrow = n_time, ncol = MCMC_sample_size + 1))
  
  # Store time values
  all_distances[, 1] <- XX[-1]
  names(all_distances) <- c("t", paste0('d', 1:MCMC_sample_size))
  
  for (j in 1:MCMC_sample_size) {
    colStart <- (j * 3) - 2  # X coordinate column
    
    # Ensure valid column indices exist
    if ((colStart + 1) > ncol(all_paths_wide)) {
      stop("colStart + 1 exceeds the column limit of all_paths_wide")
    }
    
    # Extract distance and elevation for all timstemps
    x1 <- all_paths_wide[1:n_time, colStart]
    x2 <- all_paths_wide[2:(n_time+1), colStart]
    y1 <- all_paths_wide[1:n_time, colStart + 1]
    y2 <- all_paths_wide[2:(n_time+1), colStart + 1]
    
    # Ensure vectors have the same length
    if (length(x1) != length(x2) || length(y1) != length(y2)) {
      stop("Mismatch in vector lengths")
    }
    
    # Compute slopes in a **vectorized** way
    
    all_distances[, j + 1] <- sqrt((x1 - x2)^2 + (y1 - y2)^2)
  }
  
  return(all_distances)
}

# Run the distance traveled function
all_distances <- compute_distances(all_paths_wide, MCMC_sample_size)
dim(all_distances)

# Function to compute slope between points
compute_slopes <- function(all_paths_wide, MCMC_sample_size) {
  n_time <- nrow(all_paths_wide) - 1  # Number of time steps
  
  # Create a correctly indexed data frame
  all_slopes <- data.frame(matrix(NA, nrow = n_time, ncol = MCMC_sample_size + 1))
  
  # Store time values
  all_slopes[, 1] <- XX[-1]
  names(all_slopes) <- c("t", paste0('s', 1:MCMC_sample_size))
  
  for (j in 1:MCMC_sample_size) {
    colStart <- (j * 3) - 2  # X coordinate column
    
    # Ensure valid column indices exist
    if ((colStart + 1) > ncol(all_paths_wide)) {
      stop("colStart + 1 exceeds the column limit of all_paths_wide")
    }
    
    # Extract coordinates for all time steps
    
    
    d <- all_distances[1:n_time,j]
    z1 <- all_paths_wide[1:n_time, colStart+2]
    z2 <- all_paths_wide[2:(n_time+1), colStart+2]
    
    # Ensure vectors have the same length
    if (length(z1) != length(z2)) {
      stop("Mismatch in vector lengths")
    }
    
    
    # Compute Euclidean distances in a **vectorized** way
    adjacent <- d
    opposite <- z2 - z1
    all_slopes[, j + 1] <- atan(opposite / adjacent)
  }
  
  return(all_slopes)
}

all_slopes <- compute_slopes(all_paths_wide, MCMC_sample_size)


####### NEXT THING TO DO IS CALCULATE ENERGY  ###########

# Define the ARC model function
arc_model <- function(m, theta, d) {
  8 * m^0.66 + 100 * (1 + sin(2 * theta - 74)) * m^0.88 * d
}

